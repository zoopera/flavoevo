#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Spreadsheet::WriteExcel;

my %codon = (
  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  TAT => "Y", TAC => "Y", TAA => "STOP", TAG => "STOP",
  TGT => "C", TGC => "C", TGA => "STOP", TGG => "W",
  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G");

# parse genome annotation file
my $gbff = "00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.gbf"; # use the updated new reference genome
my %ginfo = %{ginfo_extractor($gbff)};

# parse snps
my $snp_parsed_input = "SNP.parsed";
my %smp2stats = ();
open(FH_IN,"<",$snp_parsed_input) or die "Can't open $snp_parsed_input: $!";
while(my $line = <FH_IN>)
{
  chomp($line);
  my ($chr,$pos,$ref,$info,$mut_map) = split(/\t/,$line);
  my ($focal_index,$focal_max) = $mut_map =~ /([0-9]+)\|([ACGT]+)\|A\(/;
  my $type = (exists $ginfo{'feat'}{$chr}{$pos})? $ginfo{'feat'}{$chr}{$pos}{'type'}:"IGR"; # IGR = intergenic region
  if($type eq "CDS")
  {
    my $feat = (exists $ginfo{'feat'}{$chr}{$pos})? $ginfo{'feat'}{$chr}{$pos}{'info'}:"";
    my ($feat_tag,$feat_pos,$feat_nt) = split(/\|/,$feat);
    my ($feat_strand) = $feat_pos =~ /:([+-]+)/;

    my $nt_pre = $feat_nt;
    my $nt_pos = "";
    # reverse to forward
    if($feat_strand eq "-")
    {
      $nt_pre = uc(reverse($nt_pre));
      $nt_pre =~ tr/[ATGCN]/[TACGN]/;
    }
    # processing
    if($feat_pos =~ /join/)
    {
      my @subloc = $feat_pos =~ /([0-9]+\.\.[0-9]+)/g;
      my $nt_point = 0;
      for(my $i=0; $i<=$#subloc; $i++)
      {
        my ($subloc_start,$subloc_end) = $subloc[$i] =~ /([0-9]+)\.\.([0-9]+)/;
        if($pos>=$subloc_start && $pos<=$subloc_end)
        {
          # mut
          $nt_pos .= $focal_max.substr($nt_pre,$nt_point+1) if $pos==$subloc_start;
          $nt_pos .= substr($nt_pre,$nt_point,$pos-$subloc_start).$focal_max if $pos==$subloc_end;
          $nt_pos .= substr($nt_pre,$nt_point,$pos-$subloc_start).$focal_max.substr($nt_pre,$nt_point+$pos-$subloc_start+1,$subloc_end-$pos) if $pos>$subloc_start && $pos<$subloc_end;
        }else
        {
          $nt_pos .= substr($nt_pre,$nt_point,$subloc_end-$subloc_start+1);
        }
        $nt_point += $subloc_end-$subloc_start+1;
      }
    }else
    {
      my ($feat_start,$feat_end) = $feat_pos =~ /:[+-]+([0-9]+)\.\.([0-9]+)/;
      # mut
      $nt_pos = $focal_max.substr($nt_pre,1) if $pos==$feat_start;
      $nt_pos = substr($nt_pre,0,$pos-$feat_start).$focal_max if $pos==$feat_end;
      $nt_pos = substr($nt_pre,0,$pos-$feat_start).$focal_max.substr($nt_pre,$pos-$feat_start+1) if $pos>$feat_start && $pos<$feat_end;
    }
    # forward to reverse
    if($feat_strand eq "-")
    {
      $nt_pre = uc(reverse($nt_pre));
      $nt_pos = uc(reverse($nt_pos));
      $nt_pre =~ tr/[ATGCN]/[TACGN]/;
      $nt_pos =~ tr/[ATGCN]/[TACGN]/;
    }
    # trans
    my @nt_pre = split(//,$nt_pre);
    my @nt_pos = split(//,$nt_pos);
    my $codon_pre = "";
    my $codon_pos = "";
    my $aa_pre = "";
    my $aa_pos = "";
    for(my $b=0; $b<=$#nt_pre; $b+=3)
    {
      $codon_pre = join("",@nt_pre[$b..$b+2]);
      $codon_pos = join("",@nt_pos[$b..$b+2]);
      if($codon_pre ne $codon_pos) # the codon with the mut
      {
        $aa_pre = $codon{$codon_pre};
        $aa_pos = $codon{$codon_pos};
        last;
      }
    }
    # syn or non?
    my $sub_type = ($aa_pre eq $aa_pos)? "syn":"non";
    my $flg_tmp = ($mut_map =~ /^\*/)? "*":"";
    my ($map_tmp) = $mut_map =~ /^\*?[0-9]+\|(\S+)$/;
    $smp2stats{$focal_index}{$chr}{$pos} = "$ref>$focal_max\t$type\t$feat_tag ($feat_pos)\t$codon_pre>$codon_pos\t$aa_pre>$aa_pos\t$sub_type\t$flg_tmp$map_tmp";
  }elsif($type =~ /RNA/)
  {
    my $feat = (exists $ginfo{'feat'}{$chr}{$pos})? $ginfo{'feat'}{$chr}{$pos}{'info'}:"";
    my ($feat_tag,$feat_pos,$feat_nt) = split(/\|/,$feat);
    my $flg_tmp = ($mut_map =~ /^\*/)? "*":"";
    my ($map_tmp) = $mut_map =~ /^\*?[0-9]+\|(\S+)$/;
    $smp2stats{$focal_index}{$chr}{$pos} = "$ref>$focal_max\t$type\t$feat_tag ($feat_pos)\t\t\t\t$flg_tmp$map_tmp"; 
  }else
  {
    my $flg_tmp = ($mut_map =~ /^\*/)? "*":"";
    my ($map_tmp) = $mut_map =~ /^\*?[0-9]+\|(\S+)$/;
    $smp2stats{$focal_index}{$chr}{$pos} = "$ref>$focal_max\t$type\t\t\t\t\t$flg_tmp$map_tmp";
  }
}close(FH_IN) or die "Can't open $snp_parsed_input: $!";

# ===================================================================================================================================================
# output
my $workbook = Spreadsheet::WriteExcel->new("Tab2.xls");
my $worksheet = $workbook->add_worksheet();
my $format = $workbook->add_format();
$format->set_bold();
$format->set_align('vcenter');
$format->set_align('center');
my $row = 0;
my @row = ("Line","Chr","Position","Mut","Type","Gene Position","Codon Change","Amino Acid Change","Syn/Nonsyn","Read Coverage");
$worksheet->write_row($row,0,\@row,$format);
$row++;
foreach my $smp (sort keys %smp2stats)
{
  foreach my $chr (sort keys %{$smp2stats{$smp}})
  {
    foreach my $pos (sort keys %{$smp2stats{$smp}{$chr}})
    {
      @row = split(/\t/,"$smp\t$chr\t$pos\t".$smp2stats{$smp}{$chr}{$pos});
      $worksheet->write_row($row,0,\@row);
      $row++;
    }
  }
}
$workbook->close();

# ===================================================================================================================================================
# subroutines
sub ginfo_extractor
{
  my ($gbff) = @_;
  my %ginfo = ();

  my $gbff_in = Bio::SeqIO->new(-file=>"$gbff",-format=>"genbank");
  while(my $chr = $gbff_in->next_seq)
  {
    $ginfo{'gnm_size'} += $chr->length;

    foreach my $feat ($chr->get_SeqFeatures)
    {
      if($feat->primary_tag =~ /(gene|CDS|RNA)/ && !$feat->has_tag('pseudo'))
      {
        my $feat_detail = feat2detail($feat);
        if($feat->location->isa('Bio::Location::SplitLocationI'))
        {
          for my $loc ($feat->location->sub_Location)
          {
            for(my $b=$loc->start; $b<=$loc->end; $b++)
            {
              $ginfo{'feat'}{$chr->id}{$b}{'type'} = $feat->primary_tag if (!exists $ginfo{'feat'}{$chr->id}{$b} && $feat->primary_tag eq "gene") || $feat->primary_tag =~ /(CDS|RNA)/;
              $ginfo{'feat'}{$chr->id}{$b}{'info'} = $feat_detail;
            }
          }
        }else
        {
          for(my $b=$feat->start; $b<=$feat->end; $b++)
          {
            $ginfo{'feat'}{$chr->id}{$b}{'type'} = $feat->primary_tag if (!exists $ginfo{'feat'}{$chr->id}{$b} && $feat->primary_tag eq "gene") || $feat->primary_tag =~ /(CDS|RNA)/;
            $ginfo{'feat'}{$chr->id}{$b}{'info'} = $feat_detail;
          }
        }
      }
    }
  }

  return \%ginfo;
}

sub feat2detail
{
  my ($feat) = @_;
  my $feat_detail = "";

  my ($locus_tag) = ($feat->has_tag('locus_tag'))? $feat->get_tag_values('locus_tag'):"";
  my $strand = ($feat->strand eq 1)? "+":"-";
  my $subloc = "";
  if($feat->location->isa('Bio::Location::SplitLocationI')) 
  {
    for my $loc ($feat->location->sub_Location)
    {
      if($subloc eq "")
      {
        $subloc = "join(".$loc->start."..".$loc->end;
      }else
      {
        $subloc .= ",".$loc->start."..".$loc->end if $subloc ne "";
      }
    }
    $subloc .= ")";
  }else
  {
    $subloc = $feat->start."..".$feat->end;
  }
  $feat_detail = "$locus_tag|".$feat->seq_id.":$strand$subloc|".$feat->spliced_seq->seq;

  return $feat_detail;
}


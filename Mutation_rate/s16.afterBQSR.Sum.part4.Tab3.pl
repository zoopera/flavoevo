#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Spreadsheet::WriteExcel;

# parse genome annotation file
my $gbff = "00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.gbf";
my %ginfo = %{ginfo_extractor($gbff)};

# parse snps
my $indel_parsed_input = "INDEL.parsed";
my %smp2stats = ();
open(FH_IN,"<",$indel_parsed_input) or die "Can't open $indel_parsed_input: $!";
while(my $line = <FH_IN>)
{
  chomp($line);
  my ($chr,$pos,$ref,$info,$mut_map) = split(/\t/,$line);

  my ($focal_index) = $info =~ /^([0-9]+)\|/;
  my ($indel_type,$indel_base) = $mut_map =~ /^([+-]+)([^\|]+)\|/;
  my $mut = ($indel_type eq "+")? "$ref>$ref$indel_base":"$ref$indel_base>$ref";
  my $size = ($indel_type eq "+")? length($indel_base):"-".length($indel_base);
  my $type = (exists $ginfo{'feat'}{$chr}{$pos})? $ginfo{'feat'}{$chr}{$pos}{'type'}:"IGR"; # IGR = intergenic region
  if($type ne "IGR")
  {
    my $feat = (exists $ginfo{'feat'}{$chr}{$pos})? $ginfo{'feat'}{$chr}{$pos}{'info'}:"";
    my ($feat_tag,$feat_pos,$feat_nt) = split(/\|/,$feat);
    $smp2stats{$focal_index}{$chr}{$pos} = "$mut\t$size\t$type\t$feat_tag ($feat_pos)\t$mut_map";
  }else
  {
    $smp2stats{$focal_index}{$chr}{$pos} = "$mut\t$size\t$type\t\t$mut_map";
  }
}close(FH_IN) or die "Can't open $indel_parsed_input: $!";

# ===================================================================================================================================================
# output
my $workbook = Spreadsheet::WriteExcel->new("Tab3.xls");
my $worksheet = $workbook->add_worksheet();
my $format = $workbook->add_format();
$format->set_bold();
$format->set_align('vcenter');
$format->set_align('center');
my $row = 0;
my @row = ("Line","Chr","Position","Mut","Size Change","Type","Gene Position","Read Coverage");
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
        for(my $b=$feat->start; $b<=$feat->end; $b++)
        {
          $ginfo{'feat'}{$chr->id}{$b}{'type'} = $feat->primary_tag if (!exists $ginfo{'feat'}{$chr->id}{$b} && $feat->primary_tag eq "gene") || $feat->primary_tag =~ /(CDS|RNA)/;
          $ginfo{'feat'}{$chr->id}{$b}{'info'} = $feat_detail;
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


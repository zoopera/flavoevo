#!/usr/bin/env perl

use strict;
use warnings;

use File::Find::Rule;
use Spreadsheet::ParseExcel;
use Spreadsheet::WriteExcel;
use Math::Round qw(nearest);
use Bio::SeqIO;

# parse genome annotation file
my $gbff = "00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.gbf"; # use the updated reference genome
my %ginfo = %{ginfo_extractor($gbff)};

# prepare kaks_calculator input
my $dir = "00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka";
my $prot = "$dir/Sulfitobacter_sp_EE-36.faa";
my $gene = "$dir/Sulfitobacter_sp_EE-36.ffn";
my %gene = ();

my $fasta_in = Bio::SeqIO->new(-file=>"$gene",-format=>"fasta");
while(my $seq = $fasta_in->next_seq)
{
  my $sid = $seq->id;
  $sid =~ s/.*\|//g;
  $gene{$sid} = $seq->seq;
}$fasta_in->close();

open(FH_OUT,">","nt-aln.axt") or die "Can't open nt-aln.fna: $!";
$fasta_in = Bio::SeqIO->new(-file=>"$prot",-format=>"fasta");
while(my $seq = $fasta_in->next_seq)
{
  my $sid = $seq->id;
  $sid =~ s/.*\|//g;
  print FH_OUT "$sid-$sid\n".$gene{$sid}."\n".$gene{$sid}."\n\n";
}$fasta_in->close();
close(FH_OUT) or die "Can't open nt-aln.fna: $!";

# kaks_calculator
`/home-user/software/KaKs_Calculator2.0/bin/KaKs_Calculator -i nt-aln.axt -o nt-aln.kaks -m YN`;

# genomic cds s/n-sites
my $ssite = 0;
my $nsite = 0;
open(FH_IN,"<","nt-aln.kaks") or die "Can't open nt-aln.axt: $!";
my $line = <FH_IN>; # head line
while($line = <FH_IN>)
{
  chomp($line);
  my @tmp = split(/\t/,$line);
  $ssite += $tmp[7];
  $nsite += $tmp[8];
}close(FH_IN) or die "Can't close nt-aln.axt: $!";

# parse Tab2.xls
my ($tot,$smut,$nmut,$igr,$gen,$at2gc,$at2cg,$gc2at,$gc2ta) = (0,0,0,0,0,0,0,0,0);
my $parser = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse('Tab2.xls');
die $parser->error(),"\n" if !defined $workbook;
my ($worksheet) = $workbook->worksheets();
my ($row_min,$row_max) = $worksheet->row_range();
$tot = $row_max;
for my $row ($row_min .. $row_max)
{
  $smut++ if $worksheet->get_cell($row,8) &&  $worksheet->get_cell($row,8)->value() eq "syn";
  $nmut++ if $worksheet->get_cell($row,8) &&  $worksheet->get_cell($row,8)->value() eq "non";
  $igr++ if $worksheet->get_cell($row,4) && $worksheet->get_cell($row,4)->value() eq "IGR";
  $gen++ if $worksheet->get_cell($row,4) && $worksheet->get_cell($row,4)->value() ne "IGR";
  $at2gc++ if $worksheet->get_cell($row,3) && ($worksheet->get_cell($row,3)->value() eq "A>G" || $worksheet->get_cell($row,3)->value() eq "T>C");
  $at2cg++ if $worksheet->get_cell($row,3) && ($worksheet->get_cell($row,3)->value() eq "A>C" || $worksheet->get_cell($row,3)->value() eq "T>G");
  $gc2at++ if $worksheet->get_cell($row,3) && ($worksheet->get_cell($row,3)->value() eq "G>A" || $worksheet->get_cell($row,3)->value() eq "C>T");
  $gc2ta++ if $worksheet->get_cell($row,3) && ($worksheet->get_cell($row,3)->value() eq "G>T" || $worksheet->get_cell($row,3)->value() eq "C>A");
}

# parse Tab3.xls
my ($inum,$dnum,$isize,$dsize) = (0,0,0,0);
$workbook = $parser->parse('Tab3.xls');
die $parser->error(),"\n" if !defined $workbook;
($worksheet) = $workbook->worksheets();
($row_min,$row_max) = $worksheet->row_range();
for my $row ($row_min+1 .. $row_max)
{
  $inum++ if $worksheet->get_cell(0,4)->unformatted() eq "Size Change" && $worksheet->get_cell($row,4)->value() > 0;
  $dnum++ if $worksheet->get_cell(0,4)->unformatted() eq "Size Change" && $worksheet->get_cell($row,4)->value() < 0;
  $isize += $worksheet->get_cell($row,4)->value() if $worksheet->get_cell(0,4)->unformatted() eq "Size Change" && $worksheet->get_cell($row,4)->value() > 0;
  $dsize += $worksheet->get_cell($row,4)->value() if $worksheet->get_cell(0,4)->unformatted() eq "Size Change" && $worksheet->get_cell($row,4)->value() < 0;
}

# output
$workbook = Spreadsheet::WriteExcel->new("Tab4.xls");
$worksheet = $workbook->add_worksheet('STATS');
my $format = $workbook->add_format();
$format->set_bold();
$format->set_align('vcenter');
$format->set_align('center');
my $row = 0;

# part 1
my @row = ("","Base","Observation","Expectation","x2","df","p-val");
$worksheet->write_row($row,1,\@row,$format);
$row++;
@row = ("Base Substitutions",$ginfo{'gnm_size'},$tot,"","","","");
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("Intergenic",$ginfo{'gnm_size'}-$ginfo{'gen_size'},$igr,nearest(1,$tot*($ginfo{'gnm_size'}-$ginfo{'gen_size'})/$ginfo{'gnm_size'}),"","","");
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("CDS",nearest(0.01,$ssite+$nsite),($smut+$nmut),nearest(1,$tot*($ssite+$nsite)/$ginfo{'gnm_size'}),"","","");
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("Nonsynonymous Sites",nearest(0.01,$nsite),$nmut,nearest(1,$tot*$nsite/$ginfo{'gnm_size'}),"","","");
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("Synonymous Sites",nearest(0.01,$ssite),$smut,nearest(1,$tot*$ssite/$ginfo{'gnm_size'}),"","","");
$worksheet->write_row($row,1,\@row);
$row++;

# part 2
$row++;
@row = ("GC Content",100*$ginfo{'GC'}/$ginfo{'gnm_size'});
$worksheet->write_row($row,1,\@row,$format);
$row++;
@row = ("AT",$ginfo{'AT'});
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("GC",$ginfo{'GC'});
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("AT>GC && AT>CG",$at2gc+$at2cg);
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("GC>AT && GC>TA",$gc2at+$gc2ta);
$worksheet->write_row($row,1,\@row);
$row++;

# part 3
$row++;
@row = ("","Events","Bases");
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("Insertion",$inum,$isize);
$worksheet->write_row($row,1,\@row);
$row++;
@row = ("Deletion",$dnum,$dsize);
$worksheet->write_row($row,1,\@row);
$row++;

$workbook->close();

# subroutines
sub ginfo_extractor
{
  my ($gbff) = @_;
  my %ginfo = ();

  my $gbff_in = Bio::SeqIO->new(-file=>"$gbff",-format=>"genbank");
  while(my $chr = $gbff_in->next_seq)
  {
    $ginfo{'gnm_size'} += $chr->length;
    my $gc = ($chr->seq =~ tr/[GC]+//);
    my $at = ($chr->seq =~ tr/[AT]+//);
    $ginfo{'GC'} += $gc;
    $ginfo{'AT'} += $at;

    foreach my $feat ($chr->get_SeqFeatures)
    {
      if($feat->primary_tag eq "gene" && !$feat->has_tag('pseudo'))
      {
        $ginfo{'gen_size'} += $feat->length;
      }
    }
  }

  return \%ginfo;
}

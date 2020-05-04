#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(sum);
use Spreadsheet::WriteExcel;
use Statistics::Basic qw(:all);

my $gen = 4320; # (28*7+2)*24/1.1=4320
my %smp2sites = ();
my %smp2stats = ();

my $analyzed_sites_input = "AnalyzedSites.new.tab";
open(FH_IN,"<",$analyzed_sites_input) or die "Can't open $analyzed_sites_input: $!";
chomp(my $line = <FH_IN>);
while($line = <FH_IN>)
{
  chomp($line);
  my ($smp,$total,$a_cnt,$c_cnt,$g_cnt,$t_cnt,$perc_cov) = split(/\t/,$line);
  $smp2sites{$smp}{'total'} = $total;
  $smp2sites{$smp}{'A'} = $a_cnt;
  $smp2sites{$smp}{'C'} = $c_cnt;
  $smp2sites{$smp}{'G'} = $g_cnt;
  $smp2sites{$smp}{'T'} = $t_cnt;
}close(FH_IN) or die "Can't close $analyzed_sites_input: $!";

# read in snp sites
my $snp_parsed_input = "SNP.parsed";
open(FH_IN,"<",$snp_parsed_input) or die "Can't open $snp_parsed_input: $!";
while(my $line = <FH_IN>)
{
  chomp($line);
  my ($chr,$pos,$ref,$info,$mut) = split(/\t/,$line);
  my ($smp,$alt) = $mut =~ /([0-9+]+)\|([ACGT]+)\|A\(/;
  $smp2stats{$smp}{'Ts'}{'AT>GC'}++ if ($ref eq "A" && $alt eq "G") || ($ref eq "T" && $alt eq "C");
  $smp2stats{$smp}{'Ts'}{'GC>AT'}++ if ($ref eq "G" && $alt eq "A") || ($ref eq "C" && $alt eq "T");
  $smp2stats{$smp}{'Tv'}{'AT>TA'}++ if ($ref eq "A" && $alt eq "T") || ($ref eq "T" && $alt eq "A");
  $smp2stats{$smp}{'Tv'}{'GC>TA'}++ if ($ref eq "G" && $alt eq "T") || ($ref eq "C" && $alt eq "A");
  $smp2stats{$smp}{'Tv'}{'AT>CG'}++ if ($ref eq "A" && $alt eq "C") || ($ref eq "T" && $alt eq "G");
  $smp2stats{$smp}{'Tv'}{'GC>CG'}++ if ($ref eq "G" && $alt eq "C") || ($ref eq "C" && $alt eq "G");
}close(FH_IN) or die "Can't close $snp_parsed_input: $!";

# read in indel sites
my $indel_parsed_input = "INDEL.parsed";
open(FH_IN,"<",$indel_parsed_input) or die "Can't open $indel_parsed_input: $!";
while(my $line = <FH_IN>)
{
  my ($chr,$pos,$ref,$info,$mut_sig) = split(/\t/,$line);
  my ($smp) = $info =~ /^([0-9]+)\|/;
  if($mut_sig =~ /^([+-]+)[ACGT]+/)
  {
    $smp2stats{$smp}{'insertion'}++ if $1 eq "+";
    $smp2stats{$smp}{'deletion'}++ if $1 eq "-";
  }
}close(FH_IN) or die "Can't close $indel_parsed_input: $!";

# output
my $workbook = Spreadsheet::WriteExcel->new("Tab1.xls");
my $worksheet = $workbook->add_worksheet();
my $format_merge = $workbook->add_format();
$format_merge->set_bold();
$format_merge->set_align('vcenter');
$format_merge->set_align('center');
my $format = $workbook->add_format();
$format->set_bold();
$format->set_align('vcenter');
$format->set_align('center');

$worksheet->merge_range(0,0,3,0,"MA Line",$format_merge);
$worksheet->merge_range(0,1,1,6,"Substitutions",$format_merge);
$worksheet->merge_range(2,1,2,2,"Transitions",$format_merge);
$worksheet->merge_range(2,3,2,6,"Transversions",$format_merge);
$worksheet->write(3,1,"AT>GC",$format);
$worksheet->write(3,2,"GC>AT",$format);
$worksheet->write(3,3,"AT>TA",$format);
$worksheet->write(3,4,"GC>TA",$format);
$worksheet->write(3,5,"AT>CG",$format);
$worksheet->write(3,6,"GC>CG",$format);
$worksheet->merge_range(0,7,1,8,"Indels",$format_merge);
$worksheet->merge_range(2,7,3,7,"Ins.",$format_merge);
$worksheet->merge_range(2,8,3,8,"Del.",$format_merge);
$worksheet->merge_range(0,9,3,9,"Ts/Tv",$format_merge);
$worksheet->merge_range(0,10,3,10,"A Sites (bp)",$format_merge);
$worksheet->merge_range(0,11,3,11,"C Sites (bp)",$format_merge);
$worksheet->merge_range(0,12,3,12,"G Sites (bp)",$format_merge);
$worksheet->merge_range(0,13,3,13,"T Sites (bp)",$format_merge);
$worksheet->merge_range(0,14,3,14,"Sites (bp)",$format_merge);
$worksheet->merge_range(0,15,3,15,"Gen.",$format_merge);
$worksheet->merge_range(0,16,3,16,"Base-sub Rate (x10^-10)/site/gen.",$format_merge);

my ($row,$col) = (4,0);
my @snp_rate = ();
foreach my $smp (sort keys %smp2sites)
{
  my $at2gc = (exists $smp2stats{$smp}{'Ts'}{'AT>GC'})? $smp2stats{$smp}{'Ts'}{'AT>GC'}:0;
  my $gc2at = (exists $smp2stats{$smp}{'Ts'}{'GC>AT'})? $smp2stats{$smp}{'Ts'}{'GC>AT'}:0;
  my $at2ta = (exists $smp2stats{$smp}{'Tv'}{'AT>TA'})? $smp2stats{$smp}{'Tv'}{'AT>TA'}:0;
  my $gc2ta = (exists $smp2stats{$smp}{'Tv'}{'GC>TA'})? $smp2stats{$smp}{'Tv'}{'GC>TA'}:0;
  my $at2cg = (exists $smp2stats{$smp}{'Tv'}{'AT>CG'})? $smp2stats{$smp}{'Tv'}{'AT>CG'}:0;
  my $gc2cg = (exists $smp2stats{$smp}{'Tv'}{'GC>CG'})? $smp2stats{$smp}{'Tv'}{'GC>CG'}:0;
  my $ins = (exists $smp2stats{$smp}{'insertion'})? $smp2stats{$smp}{'insertion'}:0;
  my $del = (exists $smp2stats{$smp}{'deletion'})? $smp2stats{$smp}{'deletion'}:0;
  my $ts = $at2gc+$gc2at;
  my $tv = $at2ta+$gc2ta+$at2cg+$gc2cg;
  my $ts2tv = ($tv>0)? $ts/$tv:"-";
  my $asite = (exists $smp2sites{$smp}{'A'})? $smp2sites{$smp}{'A'}:0;
  my $csite = (exists $smp2sites{$smp}{'C'})? $smp2sites{$smp}{'C'}:0;
  my $gsite = (exists $smp2sites{$smp}{'G'})? $smp2sites{$smp}{'G'}:0;
  my $tsite = (exists $smp2sites{$smp}{'T'})? $smp2sites{$smp}{'T'}:0;
  my $site = (exists $smp2sites{$smp}{'total'})? $smp2sites{$smp}{'total'}:0;
  my $snp_rate = 10**10 * ($ts+$tv)/$site/$gen;
  push @snp_rate,$snp_rate;

  my @row = ($smp,$at2gc,$gc2at,$at2ta,$gc2ta,$at2cg,$gc2cg,$ins,$del,$ts2tv,$asite,$csite,$gsite,$tsite,$site,$gen,$snp_rate);
  $worksheet->write_row($row,$col,\@row);
  $row++;
}
$worksheet->write($row,15,"Mean",$format);
$worksheet->write($row,16,mean(@snp_rate),$format);
$row++;
$worksheet->write($row,15,"Std",$format);
$worksheet->write($row,16,stddev(@snp_rate),$format);
$row++;
$worksheet->write($row,15,"Sem",$format);
$worksheet->write($row,16,stddev(@snp_rate)/sqrt(scalar(@snp_rate)),$format);

$workbook->close();

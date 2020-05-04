#!/usr/bin/env perl

use strict;
use warnings;

use File::Find::Rule;

my $dir = "02_GATK";

# analyzed sites
my $out = "AnalyzedSites.new.tab";
my @tbl = File::Find::Rule->maxdepth(1)->file->name(qr/cnt.analyzed_sites.tbl/)->in($dir);
open(FH_OUT,">",$out) or die "Can't open $out: $!";
print FH_OUT "LINE_INDEX\tAll\tA\tC\tG\tT\tPercent\n";
foreach my $tbl (sort @tbl)
{
  open(FH_IN,"<",$tbl) or die "Can't open $tbl: $!";
  chomp(my @line = <FH_IN>);
  print FH_OUT "$line[$#line]\n";
  close(FH_IN) or die "Can't close $tbl: $!";
}close(FH_OUT) or die "Can't close $out: $!";

# error rate
my %err = ();
@tbl = File::Find::Rule->maxdepth(1)->file->name(qr/cnt.error_rate.tbl/)->in($dir);
foreach my $tbl (sort @tbl)
{
  open(FH_IN,"<",$tbl) or die "Can't open $tbl: $!";
  while(my $line = <FH_IN>)
  {
    chomp($line);
    my ($nt,$err,$all) = $line =~ /^([ACGT]+)\s+Err = ([0-9]+)\tAll = ([0-9]+)/;
    $err{$nt}{'err'} += $err;
    $err{$nt}{'all'} += $all;
  }close(FH_IN) or die "Can't close $tbl: $!";
}
$out = "ErrorRate.new.tab";
open(FH_OUT,">",$out) or die "Can't open $out: $!";
foreach my $nt (sort keys %err)
{
  print FH_OUT "$nt Error Rate = ". (100*$err{$nt}{'err'}/$err{$nt}{'all'}) ."%\n";
}close(FH_OUT) or die "Can't close $out: $!";

__END__
-rw-r--r-- 1 ysun admin_grp 119 May 19 15:59 02_GATK/70_CAGAGAGGATCT-AAGGAGTA_L001.clean.sort.dedup.rnd3.cnt.error_rate.tbl
-rw-r--r-- 1 ysun admin_grp  77 May 19 16:00 02_GATK/75_GCTACGCTATCT-GCGTAAGA_L001.clean.sort.dedup.rnd3.cnt.analyzed_sites.tbl
A Err = 14804	All = 24806950
C Err = 125145	All = 37083912
G Err = 127177	All = 37618694
T Err = 127177	All = 24741574


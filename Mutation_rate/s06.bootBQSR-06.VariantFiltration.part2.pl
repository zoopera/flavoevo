#!/usr/bin/env perl

use strict;
use warnings;

my $boot = 0;
my $dir = "02_GATK";
my $gatk_dir = "/home-user/software/gatk/latest";
my @type = qw(snp indel);

foreach my $type (sort @type)
{
  my $vcf_input = "$dir/Combo.rnd$boot.$type.pseudo.flt.vcf";
  my $vcf_output = "$dir/Combo.rnd$boot.$type.pseudo.flt.pass.vcf";
  open(FH_IN,"<",$vcf_input) or die "Can't open $vcf_input: $!";
  open(FH_OUT,">",$vcf_output) or die "Can't open $vcf_output: $!";
  while(my $line = <FH_IN>)
  {
    chomp($line);
    print FH_OUT "$line\n" if $line !~ /$type\_filter\tAC=/;
  }
  close(FH_IN) or die "Can't close $vcf_input: $!";
  close(FH_OUT) or die "Can't close $vcf_output: $!";

  `$gatk_dir/gatk IndexFeatureFile -F $vcf_output`;
}

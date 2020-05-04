#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::R;

my $boot = 1;
my @type = qw(snp indel);

# Main
foreach my $type (sort @type)
{
  my $vcf_input = "02_GATK/Combo.rnd$boot.$type.pseudo.vcf";
  my $stats_output = "Combo.rnd$boot.$type.pseudo.mapq.tab";
  my %stats = %{stats($vcf_input)};
  open(FH_OUT,">",$stats_output) or die "Can't open $stats_output: $!";
  foreach my $key (sort keys %stats)
  {
    print FH_OUT "$key\t",join("\t",@{$stats{$key}}),"\n";
  }close(FH_OUT) or die "Can't close $stats_output: $!";

  # R
  my $pdf_output = "Combo.rnd$boot.$type.pseudo.mapq.pdf";
  my $R = Statistics::R->new();
  $R->startR;
  $R->send(qq'dat <- read.table("$stats_output",head=F)');
  $R->send('dat <- t(dat)');
  $R->send('colnames(dat) <- dat[1,]');
  $R->send('dat <- dat[-1,]');
  $R->send('class(dat) <- "numeric"');
  $R->send('dat <- data.frame(dat)');
  $R->send(qq'pdf("$pdf_output")');
  $R->send('par(mfrow=c(2,2))');
  $R->send('hist(dat$QD,main="Histogram of QD")');
  $R->send('hist(dat$FS,main="Histogram of FS")');
  $R->send('hist(dat$MQ,main="Histogram of MQ")');
  $R->send('hist(dat$SOR,main="Histogram of SOR")');
  $R->send('dev.off()');
  $R->stopR();
}

# Subroutine
sub stats
{
  my ($vcf) = @_;
  my %stats = ();

  open(FH_IN,"<",$vcf) or die "Can't open $vcf: $!";
  while(my $line = <FH_IN>)
  {
    if($line =~ /(AC=[0-9\.]+\;AF=[0-9\.]+\;AN=[0-9\.]+\;DP=[0-9\.]+\;FS=[0-9\.]+\;MLEAC=[0-9\.]+\;MLEAF=[0-9\.]+\;MQ=[0-9\.]+\;QD=[0-9\.]+\;SOR=[0-9\.]+)/)
    {
      my @qc = split(/\;/,$1);
      foreach my $qc (sort @qc)
      {
        my ($name,$score) = split(/=/,$qc);
        push @{$stats{$name}},$score;
      }
    }
  }close(FH_IN) or die "Can't close $vcf: $!";

  return \%stats;
}

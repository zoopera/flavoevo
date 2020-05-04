#!/usr/bin/env perl

use strict;
use warnings;

use File::Find::Rule;
use Bio::SeqIO;
use Math::Round qw(nearest);

my $type = "snp";
my $output = uc($type).".parsed";

my $boot = 3;
my $dir_in = "02_GATK";
my @fvcf = File::Find::Rule->maxdepth(1)->file->name(qr/[0-9]+_.*.rnd$boot.$type.flt.vcf$/)->in($dir_in);

my %smp = ();
my %pileup = ();
my %vinfo = ();

foreach my $fvcf (sort @fvcf)
{
  print "Processing $fvcf ...\n";
  my ($smp) = $fvcf =~ /([0-9]+)\_[ACGTN]+-[ACGTN]+_/;
  $smp{$smp} = 1;

  # read in pileup file
  my $pileup = $fvcf;
  $pileup =~ s/$type.flt.vcf$/cnt/g;
  %{$pileup{$smp}} = %{pos2pileup($pileup)};

  # variants
  my $vinfo_focal = vinfo_extractor($fvcf);

  foreach my $chr (sort keys %{$vinfo_focal})
  {
    foreach my $pos (sort{$a<=>$b} keys %{$vinfo_focal->{$chr}})
    {
      foreach my $alt (sort keys %{$vinfo_focal->{$chr}{$pos}})
      {
        $vinfo{$chr}{$pos}{$alt}{$smp} = $vinfo_focal->{$chr}{$pos}{$alt};
      }
    }
  }
}

# if the same snp/indel occurs in more than 50% of the samples, we would take is as mutations occurred in the ancestral colony and delete them.
foreach my $chr (sort keys %vinfo)
{
  foreach my $pos (sort{$a<=>$b} keys %{$vinfo{$chr}})
  {
    foreach my $alt (sort keys %{$vinfo{$chr}{$pos}})
    {
      delete $vinfo{$chr}{$pos}{$alt} if scalar(keys %{$vinfo{$chr}{$pos}{$alt}}) > 0.5*scalar(keys %smp);
    }
  }
}

# output
open(FH_OUT,">",$output) or die "Can't open $output: $!";
foreach my $chr (sort keys %vinfo)
{
  foreach my $pos (sort{$a<=>$b} keys %{$vinfo{$chr}})
  {
    foreach my $alt (sort keys %{$vinfo{$chr}{$pos}})
    {
      my ($ref,$mut) = $alt =~ /(\S+) -> (\S+)/;
      foreach my $smp (sort keys %{$vinfo{$chr}{$pos}{$alt}})
      {
        print FH_OUT "$chr\t$pos\t$ref\t".$vinfo{$chr}{$pos}{$alt}{$smp}."\t$smp|$mut|".$pileup{$smp}{$chr}{$pos}."\n"; # use the mutation identified by GATK, rather than pileup stats.
      }
    }
  }
}close(FH_OUT) or die "Can't close $output: $!";


# subroutine
sub pos2pileup
{
  my ($pileup) = @_;
  my %pos2pileup = ();

  open(FH_IN,"<",$pileup) or die "Can't open $pileup: $!";
  while(my $line = <FH_IN>)
  {
    chomp($line);
    my @tmp = split(/\|/,$line);
    my ($chr,$pos,$ref,$index) = @tmp[0..3];
    $pos2pileup{$tmp[0]}{$tmp[1]} = join("|",@tmp[4..7]);
  }close(FH_IN) or die "Can't close $pileup: $!";

  return \%pos2pileup;
}

sub vinfo_extractor
{
  my ($fvcf) = @_; # vcf after VariantFiltration
  my %vinfo = ();

  my $flg = 0;
  open(FH_IN,"<",$fvcf) or die "Can't open $fvcf: $!";
  while(my $line = <FH_IN>)
  {
    chomp($line);
    $flg = 1 if $line =~ /^#CHROM/;
    if($flg == 1 && $line !~ /^#CHROM/)
    {
      my ($chr,$pos,$id,$ref,$alt,$qual,$filt,$info,$fmt,$stat) = split(/\t/,$line);
      my ($gt,$ref_ad,$alt_ad,$dp,$gq,$ref_pl,$alt_pl) = $stat =~ /([0-9]+):([0-9]+),([0-9]+):([0-9]+):([0-9]+):([0-9]+),([0-9]+)/;
      $vinfo{$chr}{$pos}{"$ref -> $alt"} = "$info||$fmt||$stat" if $filt eq "PASS" && $dp>=10 && $ref_pl>=0 && $alt_pl==0 && $alt_ad/($ref_ad+$alt_ad)>=0.8;
    }
  }close(FH_IN) or die "Can't close $fvcf: $!";

  return \%vinfo;
}



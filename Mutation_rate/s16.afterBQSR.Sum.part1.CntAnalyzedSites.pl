#!/usr/bin/env perl

use strict;
use warnings;

use File::Find::Rule;
use Math::Round qw(nearest);

my $cnt_dir = "02_GATK";
my $ref_size = 3338036;
my @cnt = File::Find::Rule->maxdepth(1)->file->name(qr/[0-9]+_.*.clean.sort.dedup.rnd3.cnt/)->in($cnt_dir);
my %err = ();

open(SITE_OUT,">","AnalyzedSites.tab") or die "Can't open AnalyzedSites.tab: $!";
print SITE_OUT "LINE_INDEX\tAll\tA\tC\tG\tT\tPercent\n";
foreach my $cnt (sort @cnt)
{
  my %analyzed_site = ();
  my $focal_index = "";

  open(CNT_IN,"<",$cnt) or die "Can't open $cnt: $!";
  while(my $line = <CNT_IN>)
  {
    chomp($line);
    my @feat = split(/\|/,$line);
    my $ref = $feat[2];
    $focal_index = $feat[3];
    my @ma = @feat[4..7];

    my %focal_snp = ();
    my $all_cnt = 0;
    for(my $i=0; $i<=$#ma; $i++)
    {
      my ($nt,$fcnt,$rcnt) = $ma[$i] =~ /([ACGT]+)\(([0-9]+):([0-9]+)\)/;
      $focal_snp{$nt}{'f'} = $fcnt;
      $focal_snp{$nt}{'r'} = $rcnt;
      $focal_snp{$nt}{'all'} = $fcnt+$rcnt;
      $all_cnt += $fcnt+$rcnt;
    }

    my ($focal_max) = sort{$focal_snp{$b}{'all'}<=>$focal_snp{$a}{'all'}} keys %focal_snp;
    my $focal_max_fcnt = (exists $focal_snp{$focal_max}{'f'})? $focal_snp{$focal_max}{'f'}:0;
    my $focal_max_rcnt = (exists $focal_snp{$focal_max}{'r'})? $focal_snp{$focal_max}{'r'}:0;
    # has focal consensus?
    if($all_cnt>0 && ($focal_max_fcnt+$focal_max_rcnt)/$all_cnt>=0.8 && $focal_max_fcnt>=3 && $focal_max_rcnt>=3)
    {
      $analyzed_site{$focal_max}++;

      # Sequencing error rates were determined using all nuclear sites that contain a consensus base call
      # that is identical to the reference nucleotide. For each individual nucleotide type (A, C, G, or T), 
      # the total number of discordant bases divided by the total number of assayed bases was used as the sequencing error rate.
      if($focal_max eq $ref)
      {
        foreach my $focal_snp (sort keys %focal_snp)
        {
          if($focal_snp ne $ref)
          {
            $err{$ref}{'err'} += $focal_snp{$focal_snp}{'f'} if exists $focal_snp{$focal_snp}{'f'};
            $err{$ref}{'err'} += $focal_snp{$focal_snp}{'r'} if exists $focal_snp{$focal_snp}{'r'};
          }
          $err{$ref}{'all'} += $focal_snp{$focal_snp}{'f'} if exists $focal_snp{$focal_snp}{'f'};
          $err{$ref}{'all'} += $focal_snp{$focal_snp}{'r'} if exists $focal_snp{$focal_snp}{'r'};
        }
      }
    }
  }close(CNT_IN) or die "Can't close $cnt: $!";

  my $Acnt = (exists $analyzed_site{'A'})? $analyzed_site{'A'}:0;
  my $Ccnt = (exists $analyzed_site{'C'})? $analyzed_site{'C'}:0;
  my $Gcnt = (exists $analyzed_site{'G'})? $analyzed_site{'G'}:0;
  my $Tcnt = (exists $analyzed_site{'T'})? $analyzed_site{'T'}:0;
  print SITE_OUT "$focal_index\t". ($Acnt+$Ccnt+$Gcnt+$Tcnt). "\t$Acnt\t$Ccnt\t$Gcnt\t$Tcnt\t". nearest(0.001,100*($Acnt+$Ccnt+$Gcnt+$Tcnt)/$ref_size) ."\n";
}close(SITE_OUT) or die "Can't close AnalyzedSites.tab: $!";

# error rate
my $ErrA = (exists $err{'A'}{'err'})? $err{'A'}{'err'}:0;
my $ErrC = (exists $err{'C'}{'err'})? $err{'C'}{'err'}:0;
my $ErrG = (exists $err{'G'}{'err'})? $err{'G'}{'err'}:0;
my $ErrT = (exists $err{'T'}{'err'})? $err{'T'}{'err'}:0;
my $AllA = (exists $err{'A'}{'all'})? $err{'A'}{'all'}:0;
my $AllC = (exists $err{'C'}{'all'})? $err{'C'}{'all'}:0;
my $AllG = (exists $err{'G'}{'all'})? $err{'G'}{'all'}:0;
my $AllT = (exists $err{'T'}{'all'})? $err{'T'}{'all'}:0;;
open(ERR_OUT,">","ErrorRate.tab") or die "Can't open ErrorRate.tab: $!";
print ERR_OUT "ErrA = $ErrA\t". nearest(0.0001,100*$ErrA/$AllA) ."%\n";
print ERR_OUT "ErrC = $ErrC\t". nearest(0.0001,100*$ErrC/$AllC) ."%\n";
print ERR_OUT "ErrG = $ErrG\t". nearest(0.0001,100*$ErrG/$AllG) ."%\n";
print ERR_OUT "ErrT = $ErrG\t". nearest(0.0001,100*$ErrT/$AllT) ."%\n";
close(ERR_OUT) or die "Can't close ErrorRate.tab: $!";

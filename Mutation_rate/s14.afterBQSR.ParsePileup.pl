#!/usr/bin/env perl

use strict;
use warnings;

my $input = $ARGV[0]; # samtools pileup
my ($line_index) = $input =~ /([0-9]+)_[ACGTN]+-[ACGTN]+_L00/;
my $output = $input;
$output =~ s/pileup$/cnt/g;

open(FH_IN,"<",$input) or die "Can't open $input: $!";
open(FH_OUT,">",$output) or die "Can't open $output: $!";
while(my $line = <FH_IN>)
{
  chomp($line);
  my ($chr,$pos,$ref,$depth,$map,$qual) = split(/\t/,$line);

  my %nt = ();
  my @map = split(//,$map);
  for(my $i=0; $i<=$#map; $i++)
  {
    # indel
    if($map[$i] eq "+" || $map[$i] eq "-")
    {
      my $flg = $map[$i];
      $i++;
      my $len = "";
      while($i<=$#map && $map[$i] =~ /[0-9]+/)
      {
        $len .= $map[$i];
        $i++;
      }
      my $indel = join("",@map[$i..$i+$len-1]);
      $nt{'indel'}{$flg.$indel}{'f'}++ if $indel =~ /[ACGT]+/;
      $nt{'indel'}{$flg.uc($indel)}{'r'}++ if $indel =~ /[acgt]+/;
      $i += $len-1; # otherwise, you will skip one base
    }elsif($map[$i] eq "^") # start of a read
    {
      $i++; # skip the character after "^"; mapping quality
    }elsif($map[$i] =~ /[ACGTacgt\.\,]+/)
    {
      $nt{'snp'}{$ref}{'f'}++ if $map[$i] eq ".";
      $nt{'snp'}{$ref}{'r'}++ if $map[$i] eq ",";
      $nt{'snp'}{uc($map[$i])}{'f'}++ if $map[$i] =~ /[ACGT]+/;
      $nt{'snp'}{uc($map[$i])}{'r'}++ if $map[$i] =~ /[acgt]+/; 
    }
  }

  # output
  print FH_OUT "$chr|$pos|$ref|$line_index|";
  # snp
  foreach my $nt (qw(A C G T))
  {
    my $fcnt = (exists $nt{'snp'}{$nt}{'f'})? $nt{'snp'}{$nt}{'f'}:0;
    my $rcnt = (exists $nt{'snp'}{$nt}{'r'})? $nt{'snp'}{$nt}{'r'}:0;
    print FH_OUT "$nt($fcnt:$rcnt)";
    print FH_OUT "|" if $nt =~ /[ACG]+/;
  }
  # indel
  foreach my $nt (sort keys %{$nt{'indel'}})
  {
    my $fcnt = (exists $nt{'indel'}{$nt}{'f'})? $nt{'indel'}{$nt}{'f'}:0;
    my $rcnt = (exists $nt{'indel'}{$nt}{'r'})? $nt{'indel'}{$nt}{'r'}:0;
    print FH_OUT "|$nt($fcnt:$rcnt)";
  }
  print FH_OUT "\n";
}
close(FH_IN) or die "Can't close $input: $!";
close(FH_OUT) or die "Can't close $output: $!";

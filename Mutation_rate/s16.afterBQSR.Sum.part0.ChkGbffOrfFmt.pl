#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

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
  GGT => "G", GGC => "G", GGA => "G", GGG => "G"
);

my $gbff = "00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.gbf";
my $gbff_in = Bio::SeqIO->new(-file=>"$gbff",-format=>"genbank");

foreach my $chr ($gbff_in->next_seq)
{
  my @cds = grep {$_->primary_tag eq 'CDS' && !$_->has_tag('pseudo')} $chr->get_SeqFeatures;

  foreach my $cds (sort @cds)
  {
    my $start = $cds->start;
    my $end = $cds->end;
    my $strand = $cds->strand;
    my $nt = $cds->spliced_seq->seq;
    my ($locus_tag) = $cds->get_tag_values('locus_tag');
    my ($aa) = $cds->get_tag_values('translation');

    if(length($nt)%3 == 0)
    {
      my $trans = "";
      my @nt = split(//,$nt);
      for(my $i=0; $i<=$#nt; $i+=3)
      {
        $trans .= $codon{join("",@nt[$i..($i+2)])} if $codon{join("",@nt[$i..($i+2)])} ne "STOP";
      }

      print "Warning: $locus_tag AA don't fit!\n$aa\n$trans\n" if substr($aa,1) ne substr($trans,1);
    }else
    {
      print "Warning: $locus_tag len=".length($nt)." %3=".(length($nt)%3)."\n";
    }
  }
}

$gbff_in->close();

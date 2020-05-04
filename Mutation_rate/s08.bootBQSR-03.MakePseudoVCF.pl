#!/usr/bin/env perl

use strict;
use warnings;

use File::Find::Rule;

my $dir = "02_GATK";
my $boot = 2;
my @vcf = File::Find::Rule->maxdepth(1)->file->name(qr/rnd$boot.vcf$/)->in($dir);
foreach my $vcf (sort @vcf)
{
	my $pseudo = $vcf;
	$pseudo =~ s/vcf$/pseudo.vcf/g;
	open(FH_IN,"<",$vcf) or die "Can't open $vcf: $!";
	open(FH_OUT,">",$pseudo) or die "Can't open $pseudo: $!";
	while(my $line = <FH_IN>)
	{
		chomp($line);
		if($line =~ /#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	smp/)
		{
			my @tmp = split(/\t/,$line);
			print FH_OUT join("\t",@tmp[0..$#tmp-1]),"\tsmp\n";
		}else
		{
			print FH_OUT "$line\n";
		}
	}
	close(FH_IN) or die "Can't close $vcf: $!";
	close(FH_OUT) or die "Can't close $pseudo: $!";
}


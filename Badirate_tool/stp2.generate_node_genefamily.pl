#!/usr/bin/perl
use warnings;

open IN,"<","badirate_CSP.out";
my @file = <IN>;
close IN;

my %hash=();
my $num=();
for (my $i=0;$i<=$#file;$i++)
{
	if ($file[$i]=~/Family Size Tree/)
	{
		for (my $j=1; $j<=$#file; $j++)
		{
			if ($file[$i+$j]=~/Total Ancestral Size/)
			{
				last;
			}
			else
			{
				$line = $file[$i+$j];
				$line =~s/;/:/;
				my ($genefamily, $nwk) = $line =~/^\t\t(\S+)\s(\(.*$)/;
				$hash{$genefamily}=1;
				while ($nwk=~/(\)\d+:)/g)
				{
					$num = $1;
					$num =~s/\)//g;
					$num =~s/://g;
					push @{$genefamily}, $num;
				}
			}
		}
	}
}

open OUT, ">>","node_genefamily";
for (my $k=0;$k<=106;$k++)
{
	my $id = $k+1;
	print OUT "$id\t";
	foreach (sort keys %hash)
	{
		if (@{$_}[$k]>0)
		{
			print OUT $_."/";
		}
	}
	$num=0;
	print OUT "\n";
}

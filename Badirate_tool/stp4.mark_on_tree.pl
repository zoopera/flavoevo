#!/usr/bin/perl
use warnings;

open IN,"<","node_summary";
while (<IN>)
{
	chomp;
	my @tmp = split(/\t/,$_);
	$hash{$tmp[0]}="'+$tmp[1]/-$tmp[2]'";
}
close IN;

open IN,"<","genome.root.node.nwk";
open OUT,">","genome.root.node.mark.nwk";
$line = <IN>;
chomp $line;
while ($line=~/\)(N\d+?):/g)
{
	$node = $1;
	$line=~s/\)\Q$node\E:/\)\Q$hash{$node}\E:/g;
}
print OUT $line."\n";
close IN;
close OUT;

#!/usr/bin/perl
use warnings;

open IN,"<","node_list_manually"; ##Manually order the ancestral node 
while (<IN>)
{
	chomp;
	push @list, $_;
}
close IN;

open IN,"<","node_genefamily";
while (<IN>)
{
	chomp;
	my @tmp = split(/\t/,$_);
	$hash{$tmp[0]}=$tmp[1];
}
close IN;

open OUT,">","node_summary";
foreach my $key (sort keys %hash)
{
	if ($key!~/^107/)
	{
		for (my $i=0;$i<=$#list;$i++)
		{
			if ($list[$i] == $key)
			{
				if ($list[$i-1]>$key)
				{
					my $rev = $hash{$key};
					my $for = $hash{$list[$i-1]};
					my @list_rev = split("/",$rev);
					my @list_for = split("/",$for);
					my %hash_rev = map{$_=>1} @list_rev;
					my %hash_for = map{$_=>1} @list_for;
					my %merge_all = map {$_ => 1} @list_rev,@list_for;
					@for_only = grep {!$hash_rev{$_}} @list_for;
					@rev_only = grep {!$hash_for{$_}} @list_rev;
					my $num_for=scalar @for_only;
					my $num_rev=scalar @rev_only;
					#===Output gain and loss lists seperately===#
					#print "N$key\tloss\t@for_only\n"; #Node loss gain
					#print "N$key\tgain\t@rev_only\n"; #Node loss gain
					print OUT "N$key\t$num_rev\t$num_for\n"; #$Node gain/loss summary
					#print "N$key\t@rev_only\n";
				}
				else
				{
					for (my $j=0;$j<=$#list;$j++)
					{
						if ($list[$i-$j]>$key) ##Find the most recent ancestor
						{
							my $rev = $hash{$key};
	        	                                my $for = $hash{$list[$i-$j]};
       		        	                        my @list_rev = split("/",$rev);
                                		        my @list_for = split("/",$for);
                 		                      	my %hash_rev = map{$_=>1} @list_rev;
                               			        my %hash_for = map{$_=>1} @list_for;
		                                        my %merge_all = map {$_ => 1} @list_rev,@list_for;
        		                                @for_only = grep {!$hash_rev{$_}} @list_for;
                       		         	        @rev_only = grep {!$hash_for{$_}} @list_rev;
                               		        	my $num_for=scalar @for_only;
							my $num_rev=scalar @rev_only;
							#===Output gain and loss lists seperately===#
                               		        	#print "N$key\tloss\t@for_only\n"; #Node loss gain
		                                        #print "N$key\tgain\t@rev_only\n"; #Node loss gain
							print OUT "N$key\t$num_rev\t$num_for\n"; #Node gain\tloss summary
							#print "N$key\t@rev_only\n";
							last;
						}
					}
				}
			}
		}
	}
}
close OUT;

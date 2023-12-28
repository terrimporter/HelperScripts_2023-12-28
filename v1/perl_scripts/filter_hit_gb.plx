#!/usr/bin/perl
#April 1, 2011 by Terri Porter
#Script to filter results in hit_lineage.txt with a list of gb
#usage $perl filter_hit_gb.plx hit_lineage.txt hit_gb.list

use strict;
use warnings;

#declare var
my $i=0;
my $filter;
my $j=0;
my $list;
my $flag=0;
my $line;

#declare array
my @lineage;
my @list;

open (LINEAGE,"<",$ARGV[0]) || die ("Error cannot read from hit_lineage.txt: $!\n");
@lineage = <LINEAGE>;
close LINEAGE;

open (LIST,"<",$ARGV[1]) || die ("Error cannot read from hit_gb.list: $!\n");
@list = <LIST>;
close LIST;

open (OUT,">>","hit_lineage.filtered") || die ("Error cannot write to hit_lineage.filtered: $!\n");

while ($list[$i]) {
	$filter = $list[$i];
	chomp $filter;

	while ($lineage[$j]) {
		$line = $lineage[$j];
		chomp $line;

		if ($flag==0) {
			if ($line =~ /$filter/) {
				print OUT $line."\n";
				$flag=1;
			}
		}
		elsif ($flag==1) {
			$j++;
			next;
		}
		$j++;
	}
	if ($flag==0) {
		print "$filter\n";#gb with no lineage
	}
	$j=0;
	$flag=0;
	$i++;
}
close OUT;

#!/usr/bin/perl
#August 27, 2013 by Terri Porter
#Script to map mothurOTU designations back to centroids
#usage perl mothurOTU_centroid_map.plx rbcLF_list.txt.filtered

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $j=0;
my $OTU;
my $count=1;
my $centroid;

#declare array
my @list;
my @line;
my @OTU;

open (LIST, "<", $ARGV[0]) || die "Error cannot open list infile: $!\n";
@list = <LIST>;
close LIST;

open (OUT, ">>", "mothurOTU_readid.map") || die "Error cannto open outfile: $!\n";

while ($list[$i]) {
	$line = $list[$i];
	chomp $line;

	@line = split(/\t/, $line);

	while ($line[$j]) {
		$OTU = $line[$j];
		
		if ($j < 2) {
			$j++;
			next;
		}
		else {
			print OUT "$count\t";
			@OTU = split(/,/,$OTU);
			$centroid = $OTU[0];
			print OUT "$centroid\n";
			$count++;
			$centroid=();
			@OTU=();
		}
		$j++;
		$OTU=();
	}
	$j=0;
	$line=();
	@line=();
	$i++;
}
$i=0;
close OUT;

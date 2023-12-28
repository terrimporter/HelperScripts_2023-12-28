#!/usr/bin/perl
#Terri Porter, Aug. 22, 2010
#Script to remove empty sequences from the seqtrim outfile, then run fasta_stats to get stats after trimming
#Usage $perl remove_empty.plx < in > out

use strict;
use warnings;

#declare variables
my $line;

#declare array
my @line;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		push (@line, $line); #add to end of array
	}
	else {
		if ($line =~ /\w+/) {
			push (@line, $line);
		}
		else {
			pop (@line);
		}
	}
}
foreach (@line) {
	$line = $_;
	print $line."\n";
}

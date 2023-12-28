#!/usr/bin/perl
#Oct.13, 2010 by Terri Porter
#Script to reformat blastall.out to make the queryID start with the queryID itself instead of the clusterID (from sorting by abundance)
#usage $perl reformat_blastall.plx < blastall.out

use strict;
use warnings;

#declare variables
my $line;
my $keep;

#declare arrays
my @line;

while (<>) {
	$line = $_;
	chop $line;
	if ($line =~ /^\d+\|\*\|G\S+/) {
		@line = split (/\|/,$line);
		$keep = $line[2];
		print "$keep\n";
	}
}


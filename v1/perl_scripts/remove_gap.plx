#!/usr/bin/perl
#Terri Porter, Aug.26, 2010
#Script to remove gap characters from a fasta file
#usage $perl remove_gap.plx < in > out

use strict;
use warnings;

#declare variables
my $line;
my $newline;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		print $line."\n";
	}
	else {
		if ($line =~ /-/) {
			$line =~ s/-//g;
			print $line."\n";
		}
		else {
			print $line."\n";
		}
	}
}

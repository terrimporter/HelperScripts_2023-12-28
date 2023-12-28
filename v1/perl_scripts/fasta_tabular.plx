#!/usr/bin/perl
#Terri Porter, Aug.25, 2010
#script to convert fasta file to tabular file for sorting in excel by ID 
#usage $perl fasta_to_tabular.plx < in > out

use warnings;
use strict;

#declare variables
my $line;
my $header;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		print "$line\t";
	}
	else {
		print "$line\n";
	}
}

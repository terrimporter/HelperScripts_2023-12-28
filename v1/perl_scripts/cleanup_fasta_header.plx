#!/usr/bin/perl
#Terri Porter, Aug.31, 2010
#Script to clean up fastas headers from uclust
#usage perl cleanup_fasta_header.plx infile

use strict;
use warnings;

#declare variables
my $outfile;
my $line;
my $ID;

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
$outfile = $ARGV[0].".cleaned";
open (OUT, '>>', $outfile) || die ("Error: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line = /^>(\w{14})/;
		$ID = $1;
		print OUT ">$ID\n";
	}
	else {
		print OUT "$line\n";
	}
}


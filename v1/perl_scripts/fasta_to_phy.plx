#!/usr/bin/perl
#Terri Porter, August 19, 2010
#Script to parse fasta file with genbak header into a file that can be turned into a phylip file in mesquite, also print mapping for bpp
#Usage $perl fasta_to_table.plx infile tempfile

use strict;
use warnings;

#declare variables
my $line;
my $gi;
my $genus;
my $species;
my $details;
my $i=0; #counter

#open fasta file
open (IN, '<', $ARGV[0]) || die ("Error: Cannot open infile $!");

#open tempfile, phylip file without the header
open (TEMP, '>>', $ARGV[1]) || die ("Error: Cannot create tempfile $!");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ /gi\|(\d+)\|gb\|.+\|\s{1}(\w+)\s{1}(\w+)\s{1}(.+)/;
		$gi = $1;
		$genus = $2;
		$species = $3;
		$details = $4;
		$i++;
		print TEMP "$gi^$i\t";
	}
	else {
		print TEMP "$line\n";
	}
}
	


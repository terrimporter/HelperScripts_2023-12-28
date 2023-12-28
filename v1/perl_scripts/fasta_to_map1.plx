#!/usr/bin/perl
#Terri Porter, August 19, 2010
#Script to parse fasta file with genbak header into a map file for bpp
#Usage $perl fasta_to_table.plx infile map_tempfile

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

#create map tempfile
open (TMP, '>>', $ARGV[1]) || die ("Error: Cannot create map_tempfile $!");

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
		print TMP "$i\t$species\n";
	}
	else {
		next;
	}
}
	


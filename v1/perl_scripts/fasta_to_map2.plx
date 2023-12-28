#!/usr/bin/perl
#Terri Porter, August 20, 2010
#Script to parse map_tempfile into a map file for bpp
#Usage $perl fasta_to_map2.plx map_tempfile map

use strict;
use warnings;

#declare variables
my $line;
my $gi;
my $genus;
my $species;
my $details;
my $i=0; #counter

#open map_tempfile
open (IN, '<', $ARGV[0]) || die ("Error: Cannot open map_tempfile $!");

#create map
open (MAP, '>>', $ARGV[1]) || die ("Error: Cannot create map file $!");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ s/mazaeus/A/) {
		print MAP "$line\n";
	}
	elsif ($line =~ s/polymnia/B/) {
		print MAP "$line\n";
	}
	elsif ($line =~ s/lysimnia/C/) {
		print MAP "$line\n";
	}
	elsif ($line =~ s/menapis/D/) {
		print MAP "$line\n";
	}
}
	


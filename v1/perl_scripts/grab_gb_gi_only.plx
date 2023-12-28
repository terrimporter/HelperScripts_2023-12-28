#!/usr/bin/perl
#Oct.5, 2010 by Terri Porter
#Script to grab accession and Genus and species from original_genbank.fasta to be used to help parse SAP results
#usage $perl grab_acc_genus_species.plx < original_genbank.fasta > original_genbank.map
#modified May 5, 2011 to parse original genbank.fasta and grab a gi and gb list
#NEW usage $perl grab_gb_gi_only.plx infile

use strict;
use warnings;

#declare var
my $line;
my $gi;
my $gb;

#declare array
my @line;

open (IN, "<",$ARGV[0]) || die ("Error cannot read from genbank.fasta: $!\n");

open (OUT,">>","gi.list") || die ("Error cannot write to gi.list: $!\n");

open (OUT2,">>","gb.list") || die ("Error cannot write to gb.list: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		@line = split(/\|/,$line);
		$gi = $line[1];
		$gb = $line[3];
		print OUT $gi."\n";
		print OUT2 $gb."\n";
	}
}
close OUT;
close OUT2;
close IN;

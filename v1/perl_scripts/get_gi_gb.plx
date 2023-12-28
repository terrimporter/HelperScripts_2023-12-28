#!/usr/bin/perl
#September 13, 2010 by Terri Porter
#Script to grab the gb number from a fasta formatted file downloaded from GenBank using an entre query via fetchseqs.plx
#Usage $perl get_gi_gb.plx genbank.fasta

use strict;
use warnings;

#declare variables
my $line;
my $gi;
my $gb;

#declare arrays
my @line;

open (OUT1,">>","gi.list") || die ("Error: $!\n");
open (OUT2, ">>", "gb.list") || die ("Error: $!\n");

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		@line = split (/\|/,$line);
		$gi = $line[1];
		$gb = $line[3];
		print OUT1 "$gi\n";
		print OUT2 "$gb\n";
	}
	else {
		next;
	}
}


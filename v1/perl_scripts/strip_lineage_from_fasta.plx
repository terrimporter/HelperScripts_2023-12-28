#!/usr/bin/perl
# Teresita M. Porter, Jan. 4, 2021
# Script to strip taxonomic lineage from a RDP classifier formatted FASTA file
# USAGE perl strip_lineage_from_fasta.plx < infile > outfile

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $line;
my $acc;

# declare array
my @line;

while (<>) {
	$line = $_;
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/\t/,$line);
		$acc = $line[0];
		print $acc."\n";
	}
	else {
		print $line."\n";
	}

}

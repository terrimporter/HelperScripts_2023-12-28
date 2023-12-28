#!/usr/bin/perl
# Jan. 12/23 by Teresita M. Porter
# Script to remove gb accession, keep lineage only
# To run vsearch derep_id (dereplicates sequences with the same id)
# USAGE perl remove_acc_from_fasta.plx testNBC.fasta > testNBC.fasta2

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $lineage;
my $infile = $ARGV[0];

# declare arrays
my @in;
my @line;

# declare hashes

# read in fasta file with gb accessions in header
open (IN, "<", $infile) || die "Error, cannot open infile\n";
@in = <IN>;
close IN;

# parse the fasta file
while ($in[$i]){
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) { #header
		@line = split(/ /,$line);
		$lineage = $line[1];
		print ">".$lineage."\n";
	}
	else {
		print $line."\n";
	}

	$i++;
}
$i=0;


#!/usr/bin/perl
# Teresita M. Porter, March 8, 2023
# Skip over cellular organisms field
# Script to reformat fasta for SINTAX
# USAGE perl fasta_for_SINTAX.plx testNBC.fasta > sintax.fasta

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $accession;
my $lineage;
my $superkingdom;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;

# declare arrays
my @in;
my @line;
my @lineage;

# declare hashes

# read infile
open (IN, "<", $ARGV[0]) || die "Error, cannot open infile\n"; 
@in = <IN>;
close IN;

# reformat file for SINTAX
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) { # header
		@line = split(/ /,$line);
		$accession = $line[0];
		$accession =~ s/>//g;
		$lineage = $line[1];
		@lineage = split(/;/,$lineage);

		$superkingdom = $lineage[1]; # here just equate to domain to accomodate bacterial outgroups later on, also superkingdom not supported by SINTAX
		$kingdom = $lineage[2];
		$phylum = $lineage[3];
		$class = $lineage[4];
		$order = $lineage[5];
		$family = $lineage[6];
		$genus = $lineage[7];
		$species = $lineage[8];
		
		print ">$accession;tax=d:$superkingdom,k:$kingdom,p:$phylum,c:$class,o:$order,f:$family,g:$genus,s:$species\n";
	}
	else {
		print $line."\n";
	}
	$i++;
}
$i=0;

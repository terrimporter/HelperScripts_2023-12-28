#!/usr/bin/perl
# Teresita M. Porter, January 12, 2023
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

		$superkingdom = $lineage[0]; # here just equate to domain to accomodate bacterial outgroups later on, also superkingdom not supported by SINTAX
		$kingdom = $lineage[1];
		$phylum = $lineage[2];
		$class = $lineage[3];
		$order = $lineage[4];
		$family = $lineage[5];
		$genus = $lineage[6];
		$species = $lineage[7];
		
		print ">$accession;tax=d:$superkingdom,k:$kingdom,p:$phylum,c:$class,o:$order,f:$family,g:$genus,s:$species\n";
	}
	else {
		$line = uc $line;
		print $line."\n";
	}
	$i++;
}
$i=0;

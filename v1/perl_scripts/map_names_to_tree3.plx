#!/usr/bin/perl
# Teresita M. Porter, Dec. 14, 2019
# edited subroutine to work with lineage field created with map_fasta_to_csv.R
# Script to add names to Newick Tree
# USAGE perl map_names_to_tree3.plx arthropoda.tree F230R_matrix.csv > arthropoda_named.tre
#
# CO1 Classifier v3 family >= 0.20 (99% correct), genus >= 0.30 (99% correct), species >= 0.70 (95% correct)

use strict;
use warnings;
use Data::Dumper;

# declare var
my $i=0;
my $line;
my $name;
my $marker_otu; # molecule_filter_PCRStep_amplicon_otu
my $SampleName;
my $SampleName_marker_otu;
my $molecule;
my $bottle;
my $molecule_bottle_marker_otu;
my $tree;

# declare array
my @tree;
my @tax;
my @line;
my @tree2;

# declare hash
my %map; #key = molecule_filter_PCRStep_amplicon_otu, value = molecule_bottle_amplicon_otu

open (TREE, "<", $ARGV[0]) || die "Error cannot open tree infile: $!\n";
@tree = <TREE>;
close TREE;

open (TAX, "<", $ARGV[1]) || die "Error cannot open taxonomy csv infile: $!\n";
@tax = <TAX>;
close TAX;

# create hash that maps marker_to lineage
# using CO1v3 cutoffs for 99% accuracy at ranks, except 95% at species rank
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	if ($i == 0) { # skip header line
		$i++;
		next;
	}
	else {
		@line = split(/,/,$line);
		$marker_otu = $line[0]; #header
		$SampleName = $line[3];
		$SampleName_marker_otu = $SampleName."_".$marker_otu;
		$molecule = $line[35];
		$bottle = $line[38];
		$molecule_bottle_marker_otu = $molecule."_".$bottle."_".$marker_otu;
		$map{$SampleName_marker_otu} = $molecule_bottle_marker_otu;
		$i++;
		next;
	}
}
$i=0;

#print Dumper(%map);

while ($tree[$i]) {
	$tree = $tree[$i];
	chomp $tree;

	while ( ($SampleName_marker_otu, $molecule_bottle_marker_otu) = each %map) {
		if (exists $map{$SampleName_marker_otu}) {
			$molecule_bottle_marker_otu = $map{$SampleName_marker_otu};
			$tree =~ s/$SampleName_marker_otu:/$molecule_bottle_marker_otu:/g;
		}
	}

	print STDOUT $tree;

	$i++;

}
$i=0;

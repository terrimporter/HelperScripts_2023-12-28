#!/usr/bin/perl
# Teresita M. Porter, April 3, 2019
# Script to add names to Newick Tree
# USAGE perl map_names_to_tree.plx insecta.tree taxonomy.csv
#
# CO1 Classifier v3 family >= 0.20 (99% correct), genus >= 0.30 (99% correct), species >= 0.70 (95% correct)

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $name;
my $otu;
my $tree;

# declare array
my @tree;
my @tax;
my @line;
my @tree2;

# declare hash
my %map; #key = otu, value = name Order_Family_Genus_Species

open (TREE, "<", $ARGV[0]) || die "Error cannot open tree infile: $!\n";
@tree = <TREE>;
close TREE;

open (TAX, "<", $ARGV[1]) || die "Error cannot open taxonomy csv infile: $!\n";
@tax = <TAX>;
close TAX;

# create hash map
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	if ($i == 0) { # skip header line
		$i++;
		next;
	}
	else {
		@line = split(/,/,$line);
		($otu, $name) = parse_name(\@line);
		$map{$otu} = $name;
	}
}
$i=0;

# map names to tree

open (OUT, ">>", "insecta_named.tree") || die "Error cannot open outfile: $!\n";

while ($tree[$i]) {
	$tree = $tree[$i];
	chomp $tree;

#	@tree2 = split(/,/,$tree);
	
#	while ($tree2[$j]) {


	while ( ($otu, $name) = each %map) {
		$tree =~ s/$otu:/$name/g;
	}

	print OUT $tree."\n";

	$i++;

}
$i=0;
close OUT;

#######################################
# new subroutine to parse high confidence assignments at each taxonomic rank

sub parse_name {

my @line = @{$_[0]};
my @marker_otu;

my $marker_otu = $line[0];
@marker_otu = split(/_/, $marker_otu);
my $marker = $marker_otu[0];
my $otu = $marker_otu[1];

my $order = $line[19];
my $oBP = $line[21];
my $family = $line[22];
my $fBP = $line[24];
my $genus = $line[25];
my $gBP = $line[27];
my $species = $line[28];
my $sBP = $line[30];

my $name = $order;

if ($fBP >= 0.20) {
	$name = $name."_".$family;
}
else {
	return ($otu, $name);
}

if ($gBP >= 0.30) {
	$name = $name."_".$genus;
}
else {
	return ($otu, $name);
}

if ($sBP >=0.70) {
	$name = $name."_".$species;
}
else {
	return ($otu, $name);
}

}

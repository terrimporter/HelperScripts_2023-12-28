#!/usr/bin/perl
# Teresita M. Porter, Sept. 17, 2019
# Script to add names to Newick Tree
# USAGE perl map_names_to_tree2.plx ephemeroptera.tree ephemeroptera.csv > ephemeroptera_named.tre
#
# CO1 Classifier v3 family >= 0.20 (99% correct), genus >= 0.30 (99% correct), species >= 0.70 (95% correct)

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $name;
my $marker_otu;
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
		($marker_otu, $name) = parse_name(\@line);
		$map{$marker_otu} = $name;
		$i++;
		next;
	}
}
$i=0;


while ($tree[$i]) {
	$tree = $tree[$i];
	chomp $tree;

	while ( ($marker_otu, $name) = each %map) {
		$tree =~ s/$marker_otu:/$name:/g;
	}

	print STDOUT $tree;

	$i++;

}
$i=0;

#######################################
# new subroutine to parse high confidence assignments at each taxonomic rank

sub parse_name {

my @line = @{$_[0]};

my $marker_otu = $line[0];

my $order = $line[20];
my $oBP = $line[22];
my $family = $line[23];
my $fBP = $line[25];
my $genus = $line[26];
my $gBP = $line[28];
my $species = $line[29];
my $sBP = $line[31];

my $name = $order;

	if ($sBP >=0.70) {
		$name = $name."_".$species;
		return ($marker_otu, $name);
	}
	else {
		if ($gBP >= 0.30) {
			$name = $name."_".$genus;
			return ($marker_otu, $name);
		}
		else {
			if ($fBP >= 0.20) {
				$name = $name."_".$family;
				return($marker_otu, $name);
			}
			else {
				return ($marker_otu, $name);
			}
		}
	}

}

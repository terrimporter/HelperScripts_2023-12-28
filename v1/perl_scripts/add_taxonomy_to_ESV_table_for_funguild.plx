#!/usr/bin/perl
# Teresita M. Porter, Jan. 20/21
# Script to add qiime-like taxonomy to an ESV table for FunGuild
# USAGE perl add_taxonomy_to_ESV_table_for_funguild.plx results.csv ESV.table

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $marker_zotu;
my $zotu;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species_sh;
my $species;
my $sh;
my $species_dec;
my $species_pct;
my $taxonomy;
my $outfile = "ESV.table.txt";
my $newline;

# declare arrays
my @results;
my @line;
my @marker_zotu;
my @species_sh;
my @table;

# declare hashes
my %zotu_tax; #key = zotu; value = $taxonomy

open (IN, "<", $ARGV[0]) || die "Error cannot open results.csv: $!\n";
@results = <IN>;
close IN;

while ($results[$i]) {
	$line = $results[$i];
	chomp $line;

	@line = split(/,/,$line);
	$marker_zotu = $line[0];
	@marker_zotu = split(/_/, $marker_zotu);
	$zotu = $marker_zotu[1];

	$kingdom = $line[8];
	$phylum = $line[11];
	$class = $line[14];
	$order = $line[17];
	$family = $line[20];
	$genus = $line[23];
	$species_sh = $line[26];
	@species_sh = split(/\|/, $species_sh);
	$species = $species_sh[0];
	$sh = $species_sh[1];
	$species_dec = $line[28];
	$species_pct = $species_dec * 100;
	$species_pct = $species_pct."%";

	$taxonomy = $species_pct."|".$species."|acc|".$sh."|rep_refs|k__".$kingdom.";p__".$phylum.";c__".$class.";o__".$order.";f__".$family.";g__".$genus.";s__".$species;

#print $taxonomy."\n"; #test

	$zotu_tax{$zotu} = $taxonomy;

	$i++;

}
$i=0;

#print Dumper(\%zotu_tax);

open (IN2, "<", $ARGV[1]) || die "Error can't open ESV.table: $!\n";
@table = <IN2>;
close IN2;

open (OUT, ">>", $outfile) || die "Error can't open outfile: $!\n";

while ($table[$i]) {
	$line = $table[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$zotu = $line[0];

	if ($i==0) {
		$newline = $line."\ttaxonomy";
	}
	else {
		if (exists $zotu_tax{$zotu}) {
			$taxonomy = $zotu_tax{$zotu};
			$newline = $line."\t$taxonomy";
		}
		else {
			print "Can't find zotu $zotu in hash\n";
		}
	}

	print OUT $newline."\n";
	$i++;

}
$i=0;

close OUT;

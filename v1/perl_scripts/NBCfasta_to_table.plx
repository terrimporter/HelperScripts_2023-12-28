#!/usr/bin/perl
#March 26, 2013 by Terri Porter
#Script to create a file that can be parsed by collapse_lineage_table.plx to create a taxonomic summary of sequences in the training set for the appendix.
#usage perl NBCfasta_to_table.plx testNBC.fasta gb_org.map.uniq

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $gbLine;
my $gb;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $lineage;
my $newline;

#declare array
my @in;
my @line;
my @gbLine;
my @map;

#declare hashes
my %table;
my %map; #indexed  by gb

open (IN, "<", $ARGV[0]) || die "Error cannot open testNBC.fasta: $!\n";
@in = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Error cannot open gb_org.map.uniq: $!\n";
@map = <IN2>;
close IN2;

#hash map
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$gb = $line[0];
#	print "map:$gb\t";
	$species = $line[1];
#	print $species."\n";

	$map{$gb} = $species;

	$i++;
	$line=();
	@line=();
	$gb=();
	$species=();
}
$i=0;

#parse out lineage from fasta
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(";", $line);
		$gbLine = $line[0];
		@gbLine = split(" ", $gbLine);
		$gb = $gbLine[0];
		$gb =~ s/^>//g;
		$kingdom = $gbLine[1];
		$phylum = $line[1];
		$class = $line[2];
		$order = $line[3];
		$family = $line[4];
		$genus = $line[5];

		$newline = $kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus;

		if (exists $map{$gb}) {
			$species = $map{$gb};
			$newline = $kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\t".$gb."\t1";
			$table{$gb} = $newline;
		}
		else {
			print "Error, $gb not found in gb_org.map.uniq\n";
		}
	}

	$i++;
	$line=();
	@line=();
	$gbLine=();
	@gbLine=();
	$gb=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$newline=();
	$species=();
}
$i=0;

#print table
open (OUT, ">>", "species.table") || die "Error cannot open outfile: $!\n";

while ( ($gb,$lineage) = each (%table) ) {
	print OUT "$lineage\n";
}
close OUT;

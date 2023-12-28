#!/usr/bin/perl
# Teresita M. Porter, May 15, 2019
# Script to create taxid.parsed file from testNBC.fasta (NCBI + BOLD merged)
# USAGE perl make_taxid_parsed.plx merged.fasta

use strict;
use warnings;

# declare var
my $outfile = "taxid.parsed";
my $i=0;
my $line;
my $header;
my $id_root;
my $id;
my $root;
my $superkingdom;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $lineage;

# declare array
my @in;
my @header;
my @id_root;

# declare hash

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) { #header
		$header = $line;
		@header = split(/;/, $header);
		$id_root = shift @header;
		@id_root = split(/ /,$id_root);
		$id = shift @id_root;
		$id =~ s/^>//g; #remove greater than sign
		$root = shift @id_root;
		$superkingdom = shift @header;
		$kingdom = shift @header;
		$phylum = shift @header;
		$class = shift @header;
		$order = shift @header;
		$family = shift @header;
		$genus = shift @header;
		$species = shift @header;
#		print $species."\n"; #test

		$lineage = $root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species;
		$lineage =~ s/_/ /g; #replace underscores with spaces
		print OUT $id."\t".$lineage."\n";
		$i+=2;
		next;
	}
	else {
		$i++;
		next;
	}
}
$i=0;
close OUT;

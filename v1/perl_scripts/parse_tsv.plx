#!/usr/bin/perl
#May 9, 2012 by Terri Porter
#Download .tsv from http://v2boldsystems.org/views/datarelease.php
#Parse file into arrays and just grab Insecta with species names 
#usage perl parse_tsv.plx iBOL_phase_2.75_COI.tsv

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $header;
my $processid;
my $phylum_reg;
my $class_reg;
my $order_reg;
my $family_reg;
my $genus_reg;
my $species_reg;
my $nucraw;
my $aminoraw;
my $entry;

#declare array
my @in;
my @header;
my @entry;

#declare hash
my %phylum;
my %class;
my %order;
my %family;
my %genus;
my %species;
my %nuc;
my %amino;


open (IN,"<",$ARGV[0]) || die "Error cannot open file: $!\n";
@in = <IN>;
close IN;

open (OUT,">>","Insecta.txt") || die "Error cannot open outfile: $!\n";
print OUT "processid\tphylum_reg\tclass_reg\torder_reg\tfamily_reg\tgenus_reg\tspecies_reg\tnucraw\taminoraw\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^processid/) {
		$header = $line;
		#print "header:$header\n";#test
		@header = split(/\t/,$header);
		$processid = $header[0];
		$phylum_reg = $header[8];
		$class_reg = $header[9];
		$order_reg = $header[10];
		$family_reg = $header[11];
		$genus_reg = $header[13];
		$species_reg = $header[14];
		$nucraw = $header[30];
		$aminoraw = $header[31];

	}

	else {
		#parse data fields
		$entry = $line;
		@entry = split(/\t/,$entry);
		$processid = $entry[0];
		$phylum_reg = $entry[8];
		$class_reg = $entry[9];
		$order_reg = $entry[10];
		$family_reg = $entry[11];
		$genus_reg = $entry[13];
		$species_reg = $entry[14];
		$nucraw = $entry[30];
		$aminoraw = $entry[31];
	
		if ($class_reg =~ /Insecta/) {
			#if ($species_reg !~ /sp\.$/) {
				print OUT "$processid\t$phylum_reg\t$class_reg\t$order_reg\t$family_reg\t$genus_reg\t$species_reg\t$nucraw\t$aminoraw\n";
			#}
		}

	}
	$i++;

}
$i=0;
close OUT;

#!/usr/bin/perl
#Feb.18,2011 by Terri Porter
#Script to compare MEGAN taxon-rank report read-id, name format with name.query from NGS_processing_part1.plx
#usage $perl MEGAN_classification_check.plx megan_species.txt megan_genus.txt name.query
#modify to check SAP classificiation instead
#usage $perl SAP_classification_check.plx clones_summary_best_species.txt clone_summary_best_genera.txt name.query
#edit pattern match to grab correct read id

use strict;
use warnings;

#var
my $i=0;
my $line;
my $readid;
my $name;
my $sap_id;
my $sap_genus1;
my $sap_species;
my $sap_genus2;
my $sap_bootstrap;
my $sap_bootstrap2;
my $j=0;
my $ref_line;
my $ref_genus;
my $ref_species;
my $genus_counter=0;
my $species_counter=0;
my $test = '';
my $genus_counter2=0;
my $k=0;
my $test2 = "sp.";
my $test3 = "cf.";
my $test4 = "aff.";

#array
my @sap;
my @ref;
my @line;
my @readid;
my @name;
my @ref_line;
my @sap2;

open (SAP,"<",$ARGV[0]) || die ("Error cannot read megan species nfile: $!\n");
@sap = <SAP>;
close SAP;

open (REF,"<",$ARGV[2]) || die ("Error cannot read name.query: $!\n");
@ref = <REF>;
close REF;

while ($sap[$i]) {
	$line = $sap[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$readid = $line[0];
	$sap_genus1 = $line[1];
	$sap_species = $line[2];
	$sap_bootstrap = $line[3];
	
	@readid = split(/\|/,$readid);
	$sap_id = $readid[1];###edit here###
	#@name = split(/ /,$name);
	#$megan_genus = $name[1];
	#$megan_species = $name[2];
	#print "megan_id=$megan_id\tmegan_genus=$megan_genus\tmegan_species=$megan_species\n";#test
	
	if ($sap_species ne $test) {
	if ($sap_species ne $test2) {
	if ($sap_species ne $test3) {
	if ($sap_species ne $test4) {
	if ($sap_bootstrap >= 95) {
		
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;

			if ($ref_line =~ /$sap_id/) {
				@ref_line = split(/\t/,$ref_line);
				$ref_genus = $ref_line[1];
				$ref_species = $ref_line[2];

				if ($ref_genus eq $sap_genus1) {
					$genus_counter++;
					#print "megan_genus:$megan_genus\tref_genus:$ref_genus\n";#test
					if ($ref_species eq $sap_species) {
						$species_counter++;
						#print "megan_species:$megan_species\tref_species:$ref_species\n";#test					
					}
					else {
						print "ref: $ref_genus $ref_species\tsap: $sap_genus1 $sap_species\n";
					}
				}
			}
			$j++;
			@ref_line = ();
			
		}
		$j=0;
	}
	}
	}
	}
	}
	else {
		$i++;
		next;
	}
	$i++;
	@line = ();
	@readid = ();
	@name = ();
	$sap_genus1 ="";
	$sap_species ="";
}

open(SAP2,"<",$ARGV[1]) || ("Error cannot read megan genus file:$!\n");
@sap2 = <SAP2>;
close SAP2;

while ($sap2[$k]) {
	$line = $sap2[$k];
	chomp $line;

	@line = split(/\t/,$line);
	$readid = $line[0];
	$sap_genus2=$line[1];
	$sap_bootstrap2 = $line[2];
	@readid = split(/\|/,$readid);
	$sap_id = $readid[1];###edit here###
	#@name = split(/ /,$name);
	#$megan_genus = $name[1];
	#$megan_species = $name[2];
	#print "megan_genus from genus file: $megan_genus\n";#test
	#print "megan_species from genus file: $megan_species\n";#test
	#if ($sap_species2 eq $test) {
	if ($sap_bootstrap2 >=95) {
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;

			if ($ref_line =~ /$sap_id/) {
				@ref_line = split(/\t/,$ref_line);
				$ref_genus = $ref_line[1];
				
				if ($ref_genus eq $sap_genus2) {
					$genus_counter2++;
					#print "megan_genus2: $megan_genus\tref_genus2:$ref_genus\n";#test
				}
				else {
					print "\nref: $ref_genus\tsap: $sap_genus2\n";
				}
			}
			$j++;
			@ref_line=();
		}
		$j=0;
	}
	else {
		$k++;
		next;
	}
	$k++;
	@line = ();
	@readid = ();
	@name = ();
	$sap_genus2="";
	#$sap_species="";
}

print "\nFrom species file: correct genera = $genus_counter\tcorrect species = $species_counter\nFrom genus file: correct genera = $genus_counter2\n";

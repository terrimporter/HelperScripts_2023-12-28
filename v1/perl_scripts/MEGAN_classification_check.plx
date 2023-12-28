#!/usr/bin/perl
#Feb.18,2011 by Terri Porter
#Script to compare MEGAN taxon-rank report read-id, name format with name.query from NGS_processing_part1.plx
#usage $perl MEGAN_classification_check.plx megan_species.txt megan_genus.txt name.query

use strict;
use warnings;

#var
my $i=0;
my $line;
my $readid;
my $name;
my $megan_id;
my $megan_genus;
my $megan_species;
my $j=0;
my $ref_line;
my $ref_genus;
my $ref_species;
my $genus_counter=0;
my $species_counter=0;
my $test = '';
my $genus_counter2=0;
my $k=0;

#array
my @megan;
my @ref;
my @line;
my @readid;
my @name;
my @ref_line;
my @megan2;

open (MEGAN,"<",$ARGV[0]) || die ("Error cannot read megan species nfile: $!\n");
@megan = <MEGAN>;
close MEGAN;

open (REF,"<",$ARGV[2]) || die ("Error cannot read name.query: $!\n");
@ref = <REF>;
close REF;

while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;
	@line = split(/,/,$line);
	$readid = $line[0];
	$name = $line[1];
	@readid = split(/\|/,$readid);
	$megan_id = $readid[1];
	@name = split(/ /,$name);
	$megan_genus = $name[1];
	$megan_species = $name[2];
	#print "megan_id=$megan_id\tmegan_genus=$megan_genus\tmegan_species=$megan_species\n";#test
	
	if ($megan_species ne $test) {
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;

			if ($ref_line =~ /$megan_id/) {
				@ref_line = split(/\t/,$ref_line);
				$ref_genus = $ref_line[1];
				$ref_species = $ref_line[2];

				if ($ref_genus eq $megan_genus) {
					$genus_counter++;
					#print "megan_genus:$megan_genus\tref_genus:$ref_genus\n";#test
					if ($ref_species eq $megan_species) {
						$species_counter++;
						#print "megan_species:$megan_species\tref_species:$ref_species\n";#test					
					}
					else {
						print "ref: $ref_genus $ref_species\tmegan: $megan_genus $megan_species\n";
					}
				}
			}
			$j++;
			@ref_line = ();
			
		}
		$j=0;
	}
	else {
		$i++;
		next;
	}
	$i++;
	@line = ();
	@readid = ();
	@name = ();
	$megan_genus ="";
	$megan_species ="";
}

open(MEGAN2,"<",$ARGV[1]) || ("Error cannot read megan genus file:$!\n");
@megan2 = <MEGAN2>;
close MEGAN2;

while ($megan2[$k]) {
	$line = $megan2[$k];
	chomp $line;

	@line = split(/,/,$line);
	$readid = $line[0];
	$name=$line[1];
	@readid = split(/\|/,$readid);
	$megan_id = $readid[1];
	@name = split(/ /,$name);
	$megan_genus = $name[1];
	$megan_species = $name[2];
	#print "megan_genus from genus file: $megan_genus\n";#test
	#print "megan_species from genus file: $megan_species\n";#test
	if ($megan_species eq $test) {
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;

			if ($ref_line =~ /$megan_id/) {
				@ref_line = split(/\t/,$ref_line);
				$ref_genus = $ref_line[1];
				
				if ($ref_genus eq $megan_genus) {
					$genus_counter2++;
					#print "megan_genus2: $megan_genus\tref_genus2:$ref_genus\n";#test
				}
				else {
					print "\nref: $ref_genus\tmegan: $megan_genus\n";
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
	$megan_genus="";
	$megan_species="";
}

print "\nFrom species file: correct genera = $genus_counter\tcorrect species = $species_counter\nFrom genus file: correct genera = $genus_counter2\n";

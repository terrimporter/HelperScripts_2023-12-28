#!/usr/bin/perl
#Feb.18,2011 by Terri Porter
#Script to compare MEGAN taxon-rank report read-id, name format with name.query from NGS_processing_part1.plx
#usage $perl MEGAN_classification_check.plx megan_species.txt megan_genus.txt name.query
#March 22, 2011 tweaked to work with parsed_lineage.txt from parse_lineage.plx and check for correct classifications to the genus, family, order, and class levels
#new usage $perl MEGAN_classification_check_part2.plx megan_readid_leaf.txt parsed_lineage.txt

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
my $family_counter=0;
my $suborder_counter=0;
my $order_counter=0;
my $subclass_counter=0;
my $class_counter=0;
my $subphylum_counter=0;
my $phylum_counter=0;
my $subkingdom_counter=0;
my $kingdom_counter=0;
my $superkingdom_counter=0;
my $ref_family;
my $ref_suborder;
my $ref_order;
my $ref_class;
my $ref_subclass;
my $ref_subphylum;
my $ref_phylum;
my $ref_subkingdom;
my $ref_kingdom;
my $ref_superkingdom;
my $flag=0;

#array
my @megan;
my @ref;
my @line;
my @readid;
my @name;
my @ref_line;
my @megan2;
my @not_classified;

open (MEGAN,"<",$ARGV[0]) || die ("Error cannot read megan species nfile: $!\n");
@megan = <MEGAN>;
close MEGAN;

open (REF,"<",$ARGV[1]) || die ("Error cannot read parsed_lineage.txt: $!\n");
@ref = <REF>;
close REF;

#look for correct genus from binomial if present
while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;
	@line = split(/,/,$line);
	$readid = $line[0];
	$name = $line[1];

	if ($name =~ /(mitosporic|incertae sedis|Not assigned|sp.|cff.|aff.|group)/) {
		$name =~ s/(mitosporic|incertae sedis|Not assigned|sp.|cff.|aff.|group)//;
	}
	
	@readid = split(/\|/,$readid);
	$megan_id = $readid[0];###gi|gb###
	@name = split(/ /,$name);###split on the space because there's a space after the comma
	$megan_genus = $name[1];
	$megan_species = $name[2];

	if ($megan_species) {
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;

			if ($ref_line =~ /$megan_id/) {
				
				@ref_line = split(/\t/,$ref_line);
				$ref_genus = $ref_line[1];
				$ref_genus =~ s/^\s+//;
#				$ref_species = $ref_line[2];

				if ($ref_genus eq $megan_genus) {
					$genus_counter++;
					$flag=1;
				}
			}
			$j++;
			@ref_line = ();
			
		}
		$j=0;
	}
	else {
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;

			if ($ref_line =~ /$megan_id/) {
#				print "megan id:$megan_id\n";#test
				@ref_line = split(/\t/,$ref_line);
				print "@ref_line\n";
				$ref_genus = $ref_line[1];
				$ref_genus =~ s/^\s+//;
				$ref_family = $ref_line[2];
				$ref_family =~ s/^\s+//;
				$ref_suborder = $ref_line[3];
				$ref_suborder =~ s/^\s+//;
				$ref_order = $ref_line[4];
				$ref_order =~ s/^\s+//;
				$ref_subclass = $ref_line[5];
				$ref_subclass =~ s/^\s+//;
				$ref_class = $ref_line[6];
				$ref_class =~ s/^\s+//;
				$ref_subphylum = $ref_line[7];
				$ref_subphylum =~ s/^\s+//;
				$ref_phylum = $ref_line[8];
				$ref_phylum =~ s/^\s+//;
				$ref_subkingdom = $ref_line[9];
				$ref_subkingdom =~ s/^\s+//;
				$ref_kingdom = $ref_line[10];
				$ref_kingdom =~ s/^\s+//;
				$ref_superkingdom = $ref_line[11];
				$ref_superkingdom =~ s/^\s+//;
				
				print "ref phylum:$ref_phylum\n";#test	
				if ($ref_genus eq $megan_genus) {
					$genus_counter2++;
					$flag=1;
				}
				elsif ($ref_family eq $megan_genus) {
					$family_counter++;
					$flag=1;
				}
				elsif ($ref_suborder eq $megan_genus) {
					$suborder_counter++;
					$flag=1;
				}
				elsif ($ref_order eq $megan_genus) {
					$order_counter++;
					$flag=1;
				}
				elsif ($ref_subclass eq $megan_genus) {
					$subclass_counter++;
					$flag=1;
				}
				elsif ($ref_class eq $megan_genus) {
					$class_counter++;
					$flag=1;
				}
				elsif ($ref_subphylum eq $megan_genus) {
					$subphylum_counter++;
					$flag=1;
				}
				elsif ($ref_phylum eq $megan_genus) {
					$phylum_counter++;
					$flag=1;
				}
				elsif ($megan_genus =~ /$ref_subkingdom/) {
					$subkingdom_counter++;
					$flag=1;
				}
				elsif ($megan_genus =~ /$ref_kingdom/) {
					$kingdom_counter++;
					$flag=1;
				}
				elsif ($ref_superkingdom eq $megan_genus) {
					$superkingdom_counter++;
					$flag=1;
				}
			}
			$j++;
			@ref_line=();
		}
		$j=0;
	}
	if ($flag==0) {
		push(@not_classified,$line);
	}
	$flag=0;
	$i++;
	@line = ();
	@readid = ();
	@name = ();
	$megan_genus="";
	$megan_species="";
}

print "\ncorrect genera from binomial = $genus_counter\ncorrect genera = $genus_counter2\ncorrect families = $family_counter\ncorrect suborders = $suborder_counter\ncorrect orders = $order_counter\ncorrect subclasses = $subclass_counter\ncorrect classes = $class_counter\ncorrect subphyla = $subphylum_counter\ncorrect phyla = $phylum_counter\ncorrect subkingdom = $subkingdom_counter\ncorrect kingdom = $kingdom_counter\ncorrect superkingdom = $superkingdom_counter\n";

my $scalar = scalar(@not_classified);
print "\n$scalar sequences not correctly classified:\n";

foreach $line (@not_classified) {
	print $line."\n";
}

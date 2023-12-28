#!/usr/bin/perl
#Feb.18,2011 by Terri Porter
#Script to compare MEGAN taxon-rank report read-id, name format with name.query from NGS_processing_part1.plx
#usage $perl MEGAN_classification_check.plx megan_species.txt megan_genus.txt name.query
#modify to check SAP classificiation instead
#usage $perl SAP_classification_check.plx clones_summary_best_species.txt clone_summary_best_genera.txt name.query
#edit pattern match to grab correct read id
#March 28, 2011 tweak to compare clone_summary_best_genera.txt with parsed_lineage.txt
#Ask user what taxonomic rank needs to be compared

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
my $sap_rank;
my $rank_to_check;
my $index='nil';
my $ref_rank;
my $rank_counter=0;

#array
my @sap;
my @ref;
my @line;
my @readid;
my @name;
my @ref_line;
my @sap2;

open (SAP,"<",$ARGV[0]) || die ("Error cannot read clone_summary_best_rank.txt infile: $!\n");
@sap = <SAP>;
close SAP;

open (REF,"<",$ARGV[1]) || die ("Error cannot read parsed_lineage.txt: $!\n");
@ref = <REF>;
close REF;

print "\nPlease one enter taxonomic rank to check\nPossible options are: genus, family, suborder, order, subclass, class, subphylum, phylum, subkingdom, kingdom, superkingdom:\n";
$rank_to_check = <STDIN>;
chomp $rank_to_check;

if ($rank_to_check eq 'genus') {
	$index=1;
}
elsif ($rank_to_check eq 'family') {
	$index=2;
}
elsif ($rank_to_check eq 'suborder') {
	$index=3;
}
elsif ($rank_to_check eq 'order') {
	$index=4;
}
elsif ($rank_to_check eq 'subclass') {
	$index=5;
}
elsif ($rank_to_check eq 'class') {
	$index=6;
}
elsif ($rank_to_check eq 'subphylum') {
	$index=7;
}
elsif ($rank_to_check eq 'phylum') {
	$index=8;
}
elsif ($rank_to_check eq 'subkingdom') {
	$index=9;
}
elsif ($rank_to_check eq 'kingdom') {
	$index=10;
}
elsif ($rank_to_check eq 'superkingdom') {
	$index=11;
}
else {
	print "Not a valid rank option\n";
}
#print "rank to check: $rank_to_check\n";#test

open (OUT,">>","trash") || die ("Error cannot write to trash: $!\n");

while ($sap[$i]) {
	$line = $sap[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$readid = $line[0];
	$sap_rank = $line[1];	
	$sap_bootstrap = $line[2];
	@readid = split(/\|/,$readid);
	$sap_id = $readid[0];###edit here format gi|gb where 0=gi 1=gb###
#	print "sap_id: $sap_id\n";#test

	if ($sap_rank ne $test) {
#		print "passed test\n";#test
	if ($sap_rank ne $test2) {
#		print "passed test2\n";#test
	if ($sap_rank ne $test3) {
#		print "passed test3\n";#test
	if ($sap_rank ne $test4) {
#		print "passed test4\n";#test
	if ($sap_bootstrap >= 95) {
#		print "passed bootstrap cutoff\n";#test
		
		while ($ref[$j]) {
			$ref_line = $ref[$j];
			chomp $ref_line;
#			print "ref_line: $ref_line\n";#test
			
			if ($ref_line =~ /^$sap_id/) {
				@ref_line = split(/\t/,$ref_line);
				$ref_rank = $ref_line[$index];
				$ref_rank =~ s/^\s{1}//;
				#$ref_species = $ref_line[2];
#				print "found sap id match\n";#test
#				print "ref_rank $ref_rank\tsap_rank $sap_rank\n";#test
				if ($ref_rank eq $sap_rank) {
					$rank_counter++;
					print OUT "ref_rank $ref_rank\tsap_rank $sap_rank\n";#test
#					print "found rank\n";#test
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
	$sap_rank ="";
}

print "\nCorrect $rank_to_check = $rank_counter\n";

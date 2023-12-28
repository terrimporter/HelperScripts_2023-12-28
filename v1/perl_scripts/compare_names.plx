#!/usr/bin/perl
#Feb.11,2011 by Terri Porter
#Script to compare hit and query names files (format gi\tgenus\tspecies)
#usage $perl compare_names.plx name.query ITS.blast.parsed.all

use strict;
use warnings;

#var
my $line;
my $i=0;
my $query_acc;
my $hit_acc;
my $hit_genus;
my $hit_species;
my $exact_match_counter=0;
my $j=0;
my $line2;
my $query_genus;
my $query_species;
my $species_match_counter=0;
my $genus_match_counter=0;

#array
my @query;
my @blast;
my @line2;

open (QUERY,"<",$ARGV[0]) || die ("Error cannot open query name file: $!\n");
@query=<QUERY>;
close QUERY;

open (BLAST,"<",$ARGV[1]) || die ("Error cannot open hit name file: $!\n");
@blast=<BLAST>;
close BLAST;

while ($blast[$i]) {
	$line = $blast[$i];
	chomp $line;
	#print $line."\n";#test
	if ($line =~ /^\w+\t\w+\|\w+\.\d+\|\s+\w+\s\S+\s+/) {
		#print "found line pattern\n";#test
		$line =~ /^(\w+)\t\w+\|(\w+)\.\d+\|\s+(\w+)\s+(\S+)\s+/;
		$query_acc = $1;
		$hit_acc = $2;
		$hit_genus = $3;
		$hit_species = $4;
		#print "query acc = $query_acc\nhit acc = $hit_acc\nhit genus = $hit_genus\nhit species = $hit_species\n";#test
		#check for exact sequence match
		if ($query_acc eq $hit_acc) {
			$exact_match_counter++;
		}
	
		#check for species match, and genus match
		while ($query[$j]) {
			$line2 = $query[$j];
			chomp $line2;
			if ($line2 =~ /^$query_acc/g) {
				@line2 = split(/\t/,$line2);
				$query_genus = $line2[1];
				$query_species = $line2[2];

				if ($query_species eq $hit_species) {
					$species_match_counter++;
				}
				if ($query_genus eq $hit_genus) {
					$genus_match_counter++;
				}
			}
			$j++;
		}
	}
	$j=0;
	$i++;
}
print "\nexact matches = $exact_match_counter\nspecies matches = $species_match_counter\ngenus matches = $genus_match_counter\n";


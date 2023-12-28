#!/usr/bin/perl
#Dec. 4, 2019 screen out 'species' that are BOLD BINs (for now)
#August 10, 2016 update hard coded names and nodes file paths and db name
#March 20, 2013 edited to make a genus_species.list for reformat_list_for_entrez_taxonomy.plx
#March 18, 2013 by Terri edited to try to grab species too...
#April 4, 2012 by Terri
#Script to use a list of taxonids from MEGAN to grab the name at each rank for lineage
#usage perl taxonomy_crawl.plx taxonid.txt

use Bio::LITE::Taxonomy::NCBI;
use strict;
use warnings;

#declare var
my $taxDB;
my $i=0;
my $taxid;
my $genus;
my $species='';#new
my $words;

#declare array
my @taxids;
my @species;

$taxDB = Bio::LITE::Taxonomy::NCBI->new(		db	=>	"nt",
												names	=>	"/home/terri/taxdmp_120419/names.dmp",
												nodes	=>	"/home/terri/taxdmp_120419/nodes.dmp");


open (IN,"<",$ARGV[0]) || die "Error cannot read in taxid infile: $!|n";
@taxids = <IN>;
close IN;

open (OUT,">>","Genus_species.txt") || die "Error cannot write to taxid.parsed: $!\n";


while ($taxids[$i]) {
	$taxid = $taxids[$i];
	chomp $taxid;
	$species = $taxDB->get_term_at_level($taxid,"species");#new

	if (defined $species && length $species > 0) {
		if ($species !~ /(sp\.|nr\.|aff\.|cf\.|BOLD)/) {
			@species = split(/ /, $species);
			$words = scalar(@species);
			if ($words == 2 ) {
				$genus = $species[0];
				$species = $species[1];
				print OUT "$genus\t$species\n";
			}
		}
	}

	$i++;
	$taxid=();
	$species='';
	@species=();
	$words=();

}
$i=0;

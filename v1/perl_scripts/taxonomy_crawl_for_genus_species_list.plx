#!/usr/bin/perl
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
my $kingdom='';
my $phylum='';
my $class='';
my $order='';
my $family='';
my $genus='';
my $i=0;
my $taxid;
my $species='';#new
my $words;

#declare array
my @taxids;
my @species;

$taxDB = Bio::LITE::Taxonomy::NCBI->new(		db	=>	"nt",
												names	=>	"/home/terri/ncbi-blast-2.9.0+/db/names.dmp",
												nodes	=>	"/home/terri/ncbi-blast-2.9.0+/db/nodes.dmp");

#print "taxDB created\n";

open (IN,"<",$ARGV[0]) || die "Error cannot read in taxid infile: $!|n";
@taxids = <IN>;
close IN;

open (OUT,">>","Genus_species.txt") || die "Error cannot write to taxid.parsed: $!\n";


while ($taxids[$i]) {
	$taxid = $taxids[$i];
	chomp $taxid;
#	print "taxid: $taxid\n";

#	$kingdom = $taxDB->get_term_at_level($taxid,"kingdom");
#	$phylum = $taxDB-> get_term_at_level($taxid,"phylum");
#	$class = $taxDB-> get_term_at_level($taxid,"class");
#	$order = $taxDB-> get_term_at_level($taxid,"order");
#	$family = $taxDB->get_term_at_level($taxid,"family"); #assignment Geospora sp. H VH-22925
#	$genus = $taxDB->get_term_at_level($taxid,"genus");
	$species = $taxDB->get_term_at_level($taxid,"species");#new

#	print OUT "$taxid\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";

#	if (defined $genus && length $genus > 0) {
		if (defined $species && length $species > 0) {
			if ($species !~ /(sp\.|nr\.|aff\.|cf\.)/) {
#				print "got to here\n";
#				print "species: $species\n";
				@species = split(/ /, $species);
				$words = scalar(@species);
#				print "words: $words\n";
				if ($words == 2 ) {
					$genus = $species[0];
					$species = $species[1];
					print OUT "$genus\t$species\n";
				}
			}
#		}
	}

	$i++;
	$taxid=();
	$kingdom='';
	$phylum='';
	$class='';
	$order='';
	$family='';
	$genus='';
	$species='';
	@species=();
	$words=();

}
$i=0;

#print "tax created\n";
#print "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\n";
#print "$kingdom\t$phylum\t$class\t$order\t$family\t$genus\n";


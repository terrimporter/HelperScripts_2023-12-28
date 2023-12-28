#!/usr/bin/perl
#April 4, 2012 by Terri
#Script to use a list of taxonids from MEGAN to grab the name at each rank for lineage
#usage perl taxonomy_crawl.plx taxonid.txt

use Bio::LITE::Taxonomy::NCBI;
use strict;
use warnings;

#declare var
my $taxDB;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $i=0;
my $taxid;
my $species;

#declare array
my @taxids;

$taxDB = Bio::LITE::Taxonomy::NCBI->new(		db	=>	"NCBI",
												names	=>	"/home/terri/taxdmp/names.dmp",
												nodes	=>	"/home/terri/taxdmp/nodes.dmp");

#print "taxDB created\n";

open (IN,"<",$ARGV[0]) || die "Error cannot read in taxid infile: $!|n";
@taxids = <IN>;
close IN;

open (OUT,">>","taxid.parsed") || die "Error cannot write to taxid.parsed: $!\n";


while ($taxids[$i]) {
	$taxid = $taxids[$i];
	chomp $taxid;
#	print "taxid: $taxid\n";

	$kingdom = $taxDB->get_term_at_level($taxid,"kingdom");
	$phylum = $taxDB-> get_term_at_level($taxid,"phylum");
	$class = $taxDB-> get_term_at_level($taxid,"class");
	$order = $taxDB-> get_term_at_level($taxid,"order");
	$family = $taxDB->get_term_at_level($taxid,"family"); #assignment Geospora sp. H VH-22925
	$genus = $taxDB->get_term_at_level($taxid,"genus");
	$species = $taxDB->get_term_at_level($taxid,"species");

	print OUT "$taxid\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";

	$i++;
	$taxid=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();

}
$i=0;

#print "tax created\n";
#print "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\n";
#print "$kingdom\t$phylum\t$class\t$order\t$family\t$genus\n";


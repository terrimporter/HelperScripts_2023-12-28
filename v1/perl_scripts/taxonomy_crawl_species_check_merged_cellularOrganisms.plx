#!/usr/bin/perl
# May 7, 2019 update path to blast 2.9.0/db
#April 10, 2018 clean up script a bit; add 'cellularOrganisms' to output
#Oct. 24, 2016 updated to grab superkingdom for prokaryotes
#August 5, 2016 updated to account for MEGAN -1 (no BLAST hits) and -2 (not assigned)
#March 11, 2014 updated to process merged.dmp too
#USAGE NEW perl taxonomy_crawl_species.plx filtered_taxids.txt.uniq

#March 18, 2013 by Terri edited to try to grab species too...
#April 4, 2012 by Terri
#Script to use a list of taxonids from MEGAN to grab the name at each rank for lineage
#usage perl taxonomy_crawl.plx taxonid.txt

use Bio::LITE::Taxonomy::NCBI;
use strict;
use warnings;

#declare var
my $taxDB;
my $superkingdom;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $i=0;
my $taxid;
my $species;#new
my $merged = "/home/terri/ncbi-blast-2.9.0+/db/merged.dmp"; ### if there are taxids that are not found, then check this map file
my $old_taxid;
my $new_taxid;
my $line;

#declare array
my @taxids;
my @line;
my @merged;

#declare hash
my %merged; #key=old_taxid; value=new_taxid

$taxDB = Bio::LITE::Taxonomy::NCBI->new(		db	=>	"NCBI",
												names	=>	"/home/terri/ncbi-blast-2.9.0+/db/names.dmp",
												nodes	=>	"/home/terri/ncbi-blast-2.9.0+/db/nodes.dmp"); 
### if taxids are not found, update these files from  ncbi ftp site for taxdump.tar.gz

#process merged.dmp map file
open (MERGED, "<", $merged) || die "Error cannot open merged.dmp: $!\n";
@merged = <MERGED>;
close MERGED;

while ($merged[$i]) {
	$line = $merged[$i];
	chomp $line;

	$line =~ s/\s+//g; #remove whitespace

	if ($line =~ /^\d+\|\d+\|/) {
		$line =~ /^(\d+)\|(\d+)\|/;
		$old_taxid = $1;
		$new_taxid = $2;
		$merged{$old_taxid} = $new_taxid;
	}

	$i++;
	$line=();
	@line=();
	$old_taxid=();
	$new_taxid=();
}
$i=0;

#process existing taxid file
open (IN,"<",$ARGV[0]) || die "Error cannot read in taxid infile: $!/n";
@taxids = <IN>;
close IN;

open (OUT,">>","taxid.parsed") || die "Error cannot write to taxid.parsed: $!\n";


while ($taxids[$i]) {
	$taxid = $taxids[$i];
	chomp $taxid;

	# Relic from when I used to parse taxid lists from MEGAN
	if ($taxid == -1 ) {
		print OUT "$taxid\tNo BLAST hits\n";
	}
	
	# Relic from when I used to parse taxid lists from MEGAN
	elsif ($taxid == -2) {
		print OUT "$taxid\tNot assigned\n";
	}

	else {
		
		$superkingdom = $taxDB->get_term_at_level($taxid,"superkingdom");

		if (length $superkingdom > 0) {
			$superkingdom = $taxDB-> get_term_at_level($taxid,"superkingdom");
			$kingdom = $taxDB-> get_term_at_level($taxid,"kingdom");
			$phylum = $taxDB-> get_term_at_level($taxid,"phylum");
			$class = $taxDB-> get_term_at_level($taxid,"class");
			$order = $taxDB-> get_term_at_level($taxid,"order");
			$family = $taxDB-> get_term_at_level($taxid,"family"); #assignment Geospora sp. H VH-22925
			$genus = $taxDB-> get_term_at_level($taxid,"genus");
			$species = $taxDB-> get_term_at_level($taxid,"species");#new
			print OUT "$taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";

		}
		else {
			if (exists $merged{$taxid}) {
				$new_taxid = $merged{$taxid};
				$superkingdom = $taxDB-> get_term_at_level($new_taxid,"superkingdom");
				$kingdom = $taxDB-> get_term_at_level($new_taxid,"kingdom");
				$phylum = $taxDB-> get_term_at_level($new_taxid,"phylum");
				$class = $taxDB-> get_term_at_level($new_taxid,"class");
				$order = $taxDB-> get_term_at_level($new_taxid,"order");
				$family = $taxDB-> get_term_at_level($new_taxid,"family"); #assignment Geospora sp. H VH-22925
				$genus = $taxDB-> get_term_at_level($new_taxid,"genus");
				$species = $taxDB-> get_term_at_level($new_taxid,"species");#new
				print OUT "$new_taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
				print OUT "$taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
			}
			else {
				print "Error cannot find taxonomic information for $taxid\n";
			}
		}
	}

	$i++;
	$taxid=();
	$superkingdom=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
	$new_taxid=();

}
$i=0;

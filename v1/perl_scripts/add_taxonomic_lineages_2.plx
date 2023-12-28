#!/usr/bin/perl
# July 27, 2021 add a check for species level identification
# July 4, 2020 based on taxonomy_crawl_species_check_merged_cellularOrganisms.plx
# always keep names & nodes.dmp updated
# Usage perl add_taxonomic_lineages.plx taxonid.txt

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
my $species_parts;
my $flag;

#declare array
my @taxids;
my @line;
my @merged;
my @species;

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

		
	$superkingdom = $taxDB->get_term_at_level($taxid,"superkingdom");

	if (length $superkingdom > 0) {
		$superkingdom = $taxDB-> get_term_at_level($taxid,"superkingdom");
		$kingdom = $taxDB-> get_term_at_level($taxid,"kingdom");
		$phylum = $taxDB-> get_term_at_level($taxid,"phylum");
		$class = $taxDB-> get_term_at_level($taxid,"class");
		$order = $taxDB-> get_term_at_level($taxid,"order");
		$family = $taxDB-> get_term_at_level($taxid,"family");
		$genus = $taxDB-> get_term_at_level($taxid,"genus");
		$species = $taxDB-> get_term_at_level($taxid,"species");

		$flag = check_species_rank();
		if ($flag == 0) {	
			print OUT "$taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
		}
		else {
			print "Taxid is not identified to the species rank\n";
		}
	}
	else {
		if (exists $merged{$taxid}) {
			$new_taxid = $merged{$taxid};
			$superkingdom = $taxDB-> get_term_at_level($new_taxid,"superkingdom");
			$kingdom = $taxDB-> get_term_at_level($new_taxid,"kingdom");
			$phylum = $taxDB-> get_term_at_level($new_taxid,"phylum");
			$class = $taxDB-> get_term_at_level($new_taxid,"class");
			$order = $taxDB-> get_term_at_level($new_taxid,"order");
			$family = $taxDB-> get_term_at_level($new_taxid,"family");
			$genus = $taxDB-> get_term_at_level($new_taxid,"genus");
			$species = $taxDB-> get_term_at_level($new_taxid,"species");

			$flag = check_species_rank();
			if ($flag == 0) {
				print OUT "$new_taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
#		print OUT "$taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
			}
			else {
				print "Taxid is not identified to the species rank\n";
			}
		}
		else {
			print "Error cannot find taxonomic information for $taxid\n";
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
	$flag = 0;

}
$i=0;

#############################
sub check_species_rank {

	if ($species =~ /(sp\.|cf\.|aff\.|)/) {
		$flag = 1;
	}
	if ($species =~ /[0-9]/) {
		$flag = 1;
	}
	@species = split(/ /, $species);
	$species_parts = scalar(@species);
	if ($species_parts > 2) {
		$flag = 1;
	}

}

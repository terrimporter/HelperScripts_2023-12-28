#!/usr/bin/perl
# March 29/23 Edit to handle commas in the name_unique field
# Jan. 10, 2023 update path to taxdmp_date/, check taxdmp to keep unique name field from names.dmp to avoid having to make manual corrections when training the classifier
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

#use Bio::LITE::Taxonomy::NCBI;
use Data::Dumper;

use strict;
use warnings;

#declare var
my $taxDB;
my $superkingdom='';
my $kingdom='';
my $phylum='';
my $class='';
my $order='';
my $family='';
my $genus='';
my $species='';
my $i=0;
my $j=0;
my $taxid;
my $merged = "/mnt/usr_home/terri/taxdmp_2023-01-10/merged.dmp"; ### if there are taxids that are not found, then check this map file
my $names = "/mnt/usr_home/terri/taxdmp_2023-01-10/names.dmp"; # try to use the unique name field if present
my $nodes = "/mnt/usr_home/terri/taxdmp_2023-01-10/nodes.dmp"; # to crawl up taxonomy tree keep
my $parent_taxid; 
my $name_txt;
my $name_unique;
my $name_class;
my $old_taxid;
my $new_taxid;
my $line;
my $taxonomy;
my $node;
my $name;
my $rank;
my $hash_ref;
my $scalar_ref;
my $original_taxid;

#declare array
my @taxids;
my @line;
my @merged;
my @names;
my @nodes;
my @node;
my @arr;

#declare hash
my %merged; #key=old_taxid; value=new_taxid
my %names; #key=name_txt; value = unique_name; turn spaces into _'s
my %nodes; #key=taxid; value = name|rank|parent_taxid
my %ranks = (superkingdom => '1',
		kingdom => '1',
		phylum => '1',
		class => '1',
		order => '1',
		family => '1',
		genus => '1',
		species => '1');
my %record; #key=rank; value = name;




# hash merged.dmp map file
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

#print "Finished parsing merged...\n";

# create a new outfile to contain the curated lineages
open (OUT,">>","taxid.parsed") || die "Error cannot write to taxid.parsed: $!\n";

# hash names.dmp map file
open (NAMES, "<", $names) || die "Error cannot open names.dmp: $!\n";
@names = <NAMES>;
close NAMES;

while ($names[$i]) {
	$line = $names[$i];
	chomp $line;

	@line = split(/\|/,$line);
	$taxid = $line[0];
	$taxid =~ s/\t//g;
	$name_txt = $line[1];
	$name_txt =~ s/\t//g;
	$name_unique = $line[2];
	$name_unique =~ s/\t//g;
	$name_unique =~ s/\s+/_/g;
	$name_unique =~ s/\,/_/g; # this comma messes up the csv files!
	$name_unique =~ s/\>//g;
	$name_unique =~ s/\<//g;
	$name_class = $line[3];
	$name_class =~ s/\t//g;

	# if a unique name exists, use it
	if ($name_class eq "scientific name") {
		if (length $name_unique) {
			$names{$taxid} = $name_unique;
		}
		else {
			$names{$taxid} = $name_txt;
		}
	}

	$i++;
}
$i=0;


# hash nodes.dmp to crawl up taxonomy tree using taxids
open (NODES, "<", $nodes) || die "Error cannot open names.dmp: $!\n";
@nodes = <NODES>;
close NODES;

while ($nodes[$i]) {
	$line = $nodes[$i];
	chomp $line; # taxid | parent tax_id | rank

	@line = split(/\|/,$line);
	$taxid = $line[0];
	$taxid =~ s/\t//g;
	$parent_taxid = $line[1];
	$parent_taxid =~ s/\t//g;
	$rank = $line[2];
	$rank =~ s/\t//g;

	$nodes{$taxid} = $parent_taxid."|".$rank;

	$i++;
}
$i=0;

# process existing taxid file
open (IN,"<",$ARGV[0]) || die "Error cannot read in taxid infile: $!/n";
@taxids = <IN>;
close IN;

while ($taxids[$i]) {
	$taxid = $taxids[$i];
	chomp $taxid; # should represent a species
	$original_taxid = $taxid;

	$rank = '';

	# crawl up the taxonomy tree from species -> superkingdom, grabbing names and ranks 
	until ($rank eq 'superkingdom' || $rank eq 'NULL') { #runs block if condition is false

		if ($j < 50) { # to avoid infinite loops

			# check that I have a valid taxid, get parent taxid, get current taxid rank
			@arr = taxonomy_crawl(\%record, \$taxid);
			$hash_ref = $arr[0];
			$scalar_ref = $arr[1];
			%record = %$hash_ref; # dereference
			$taxid = $$scalar_ref; # dereference
			$j++;
		}
		else {
			print "Problem parsing taxid $taxid\n";
			$rank = 'NULL';
		}

	}
	$j=0;

#	print Dumper(\%record);

	if (exists $record{'superkingdom'}) {
		$superkingdom = $record{'superkingdom'};
	}
	if (exists $record{'kingdom'}){
		$kingdom = $record{'kingdom'};
	}
	if (exists $record{'phylum'}){
		$phylum = $record{'phylum'};
	}
	if (exists $record{'class'}){
		$class = $record{'class'};
	}
	if (exists $record{'order'}){
		$order = $record{'order'};
	}
	if (exists $record{'family'}){
		$family = $record{'family'};
	}
	if (exists $record{'genus'}){
		$genus = $record{'genus'};
	}
	if (exists $record{'species'}){
		$species = $record{'species'};
	}
	
#	print "$taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
	print OUT "$original_taxid\tcellularOrganisms\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";

	$i++;
	$taxid=();
	$superkingdom='';
	$kingdom='';
	$phylum='';
	$class='';
	$order='';
	$family='';
	$genus='';
	$species='';
	%record=();


}
$i=0;

########################

sub taxonomy_crawl {

# declare local vars
my %hash = %{$_[0]};
my $taxid = ${$_[1]};

# get parent_taxid and name (indexed by rank)
if (exists $nodes{$taxid}) {
	$node = $nodes{$taxid};
	@node = split(/\|/,$node);
	$parent_taxid = $node[0];
	$rank = $node[1];

	if (exists $ranks{$rank}){
		# recall that %names only includes unique names
		# do this to reduce having to manually edit script 
		# to accomodate the use of non-unique names across 
		# different fields of study
		$name = $names{$taxid};
		$hash{$rank} = $name;
	}
}

# these taxids may have been merged into a newer taxid
elsif (exists $merged{$taxid}) {
	$new_taxid = $merged{$taxid};
	$taxid = $new_taxid;

	if (exists $nodes{$taxid}) {
		$node = $nodes{$taxid};
		@node = split(/\|/,$node);
		$parent_taxid = $node[0];
		$rank = $node[1];

		if (exists $ranks{$rank}){
			$name = $names{$taxid};
			$hash{$rank} = $name;
		}

	}
	else { # just can't process these taxids
		print "Error cannot find taxonomic information for $taxid\n";
		$rank = 'NULL';
	}
}

# just can't process these taxids
else {
	print "Error cannot find taxonomic information for $taxid\n";
	$rank = 'NULL';
	}
#	print Dumper(\%record);
	return (\%hash, \$parent_taxid);

}

#!/usr/bin/perl
# May 10, 2019 by Teresita M. Porter
# Script to grab NCBI info for each BOLD species and create files needed to make fasta file
## filtered taxids
## id-seq
## taxid lineage from root to species with spaces (not underscores)
## id-taxid
# only work with good species names with complete lineages for now
# ignore unnamed BINS for now
# create fasta file with names like undef and taxon_undef where needed
# also create taxid.parsed
## Hard-coded path to names.dmp and nodes.dmp, adjust as needed
# Usage perl consolidate_BOLD_NCBI_taxonomy.plx BOLD_DR_species.fasta

use strict;
use warnings;

# declare var
my $namespath = "/home/terri/ncbi-blast-2.9.0+/db/names.dmp";
my $nodespath = "/home/terri/ncbi-blast-2.9.0+/db/nodes.dmp";
my $outfile = "BOLD_DR_consolidated.fasta";
my $outfile2 = "taxid.parsed"
my $i=0;
my $line;
my $taxid;
my $original_taxid;
my $taxid_list;
my $name;
my $nameclass; #only work with scientific names (ignore common)
my $parent;
my $rank;
my $species;
my $genus;
my $family;
my $order;
my $class;
my $phylum;
my $kingdom;
my $superkingdom;
my $id_root;
my $lineage;
my $count; #to get out of any recursions that go too deep
my $j;
my $seq;
my $length;

# declare array
my @names;
my @nodes;
my @line;
my @in;
my @species;

# declare hash
my %taxid_name; #key = taxid, value = name
my %name_taxid; #key = name, value = taxid
my %taxid_parent; #key = taxid, value = parent taxid
my %taxid_rank; #key = taxid, value = taxid rank

open (NAMES, "<", $namespath) || die "Cannot open names.dmp: $!\n";
@names = <NAMES>;
close NAMES;

open (NODES, "<", $nodespath) || die "Cannot open nodes.dmp: $!\n";
@nodes = <NODES>;
close NODES;

# hash names for quick lookups
while ($names[$i]) {
	$line = $names[$i];
	chomp $line;

	@line = split(/\t\|\t/,$line);
	$taxid = $line[0];
	$name = $line[1];
	$nameclass = $line[3];
	$nameclass =~ s/\t\|//g; #remove final pipe

	if ($nameclass eq "scientific name") {
		$taxid_name{$taxid} = $name; # should always be unique
		if (exists $name_taxid{$name}){ # name could have more than one taxid, do not overwrite
			$original_taxid = $name_taxid{$name};
			$taxid_list = join(";", $original_taxid, $taxid);
			$name_taxid{$name} = $taxid_list;
		}
		else {
			$name_taxid{$name} = $taxid;
		}
	}

	$i++;
}
$i=0;


# hash nodes for quick lookups
while ($nodes[$i]) {
	$line = $nodes[$i];
	chomp $line;

	@line = split(/\t\|\t/,$line);
	$taxid = $line[0];
	$parent = $line[1];
	$rank = $line[2];

	$taxid_parent{$taxid} = $parent;
	$taxid_rank{$taxid} = $rank;

	$i++;
}
$i=0;

open (IN, "<", $ARGV[0]) || die "Error can't open infile: $!\n";
@in = <IN>;
close IN;

open (FAS, ">>", $outfile) || die "Error can't open outfile: $!\n"; # fasta file

open (TAX, ">>", $outfile2) || die "Error can't open 2nd outfile: $!\n"; # taxid.parsed

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/;/,$line);
		$species = pop(@line); #grab species at end or array
		$genus = pop(@line); 
		$family = pop(@line); #consolidate BOLD family with NCBI taxonomy
		$order = pop(@line);
		$class = pop(@line);
		$phylum = pop(@line);
		$kingdom = pop(@line);
		$superkingdom = pop(@line);
		$id_root = pop(@line);

		# parse from superkingdom to species
		if (exists $name_taxid{$superkingdom}) {
			@superkingdom = split(/;/,$superkingdom);
			$length = scalar(@superkingdom);
			if ($length > 1) { #more than one taxid for this name
				#handle this problem
				print "more than one taxid for $superkingdom\n"; #deal with this on case by case basis here
			}
		}
		else {
			# BOLD superkingdom not in NCBI
			$superkingdom = "undef";
		}

		if (exists $name_taxid{$kingdom}) {
			@kingdom = split(/;/,$kingdom);
			$length = scalar(@kingdom);
			if ($length >1) { #more than one taxid for this name
				#handle this problem
				print "more than one taxid for $kingdom\n";	 # deal with this on case by case basis here
			}
		}
		else {
			$kingdom = "undef";
		}

		if (exists $name_taxid{$phylum}) {
			@phylum = split(/;/,$phylum);
			$length = scalar(@phylum);
			if ($length > 1) {
				#handle this problem
				print "more than one taxid for $phylum\n"; # deal with this on case by case basis here
			}
		}
		else {
			$phylum = "undef";
		}

		if (exists $name_taxid{$class}) {
			@class = split(/;/,$class);
			$length = scalar(@class);
			if ($length > 1) {
				print "more than one taxid for #class\n"; #deal with this on case by case basis here
			}
		}
		else {
			$class = "undef";
		}

		unless (exists $name_taxid{$order}) {
			$order = $order."_".$class;
		}

		unless (exists $name_taxid{$family}) {
			$family = $family."_".$order;
		}

		unless (exists $name_taxid{$genus}) {
			$genus = $genus."_".$family;
		}

		unless (exists $name_taxid{$species}) {
			$species = $species."_".$genus;
		}

		$lineage = $superkingdom.";".$kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
		$lineage =~ s/ /_/g; #replace spaces with underscores
		print FAS $id_root.";".$lineage."\n";
		$j = $i+1;
		$seq = $in[$j];
		chomp $seq;
		$seq = lc $seq; # make lower-case to match RDP style
		print FAS $seq."\n";
		$i+=2;
		next;
	}
	else {
		$i++;
		next;
	}
}
$i=0;
close OUT;

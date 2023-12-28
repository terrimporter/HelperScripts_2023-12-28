#!lusr/bin/perl
# Teresita M. Porter, March 10, 2023
# Script to parse through GBIF iNaturalist research grade dataset for Chordata from USA & Canada
# Check for synonyms using recent NCBI taxonomy
# Parse through latest testNBC.taxonomy file to assess coverage
# USAGE perl check_species_coverage.plx gbif_chordata_2023-03-09.csv testNBC_2.taxonomy

use strict;
use warnings;
use Data::Dumper;
use Text::CSV;

# declare vars
my $i=0;
#my $line;
my $flag=0;
my $species;
my $num;
my $names = "/mnt/usr_home/terri/taxdmp_2023-03-09/names.dmp";
my $nodes = "/mnt/usr_home/terri/taxdmp_2023-03-09/nodes.dmp";
my $taxid;
my $name_txt;
my $name_class;
my $j=0;
my $csv;
my $rank;
my $k=0;
my $l=0;

# declare arrays
my @line;
my @gbif;
my @rdp;
my @names;
my @nodes;
my @csv;

# declare hashes
my %gbif; # key = species; value = 1
my %rdp; # key = species; value = 1
my %taxids; # key = taxid; value = 1
my %scientificNames; # key = taxid; value = name_txt
my %scientificNames2; # key = species; value = taxid
my %synonyms; # key = name_txt; value = taxid

# read in names.dmp
open (NAMES, "<", $names) || die "Error can't open names.dmp file: $!\n";
@names = <NAMES>;
close NAMES;

# read in nodes.dmp
open (NODES, "<", $nodes) || die "Error can't open nodes.dmp file: $!\n";
@nodes = <NODES>;
close NODES;

# read in rdp classifier taxonomy file
open (RDP, "<", $ARGV[1]) || die "Error can't open RDP taxonomy file: $!\n";
@rdp = <RDP>;
close RDP;
$num = scalar @rdp;
#print "There are ", $num, " taxa in the RDP file\n";
$num=();

# parse gbif to get species names
# messed up because there are commas inside the fields too, ugh
$csv = Text::CSV -> new ({binary => 1, 
						auto_diag => 1, 
				allow_whitespace => 1});

my $file = $ARGV[0];
open my $fh, "<", $file or die "Cannot use CSV: " . Text::CSV->error_diag();

# skip over header line
my @headers = @{ $csv->getline($fh) };

#while ($gbif[$i]) {
while (my $line = $csv -> getline ($fh)) {
	$species = $line->[19];
	$species =~ s/ /_/;
	$gbif{$species} = 1;
	# print "GBIF species ",$species, "\n";
}


# count total number of GBIF Chordata USA & Canada species
$num = keys %gbif;

print "There are ", $num, " unique species in GBIF.\n";
$num=();

# parse the rdp species
while ($rdp[$i]) {
	my $line = $rdp[$i];
	chomp $line;

	if ($line =~ m/Chordata/ & $line =~ m/phylum$/)  {
#		print "Found Chordata\n";
#		print $line."\n";
		$flag = 1;
		$i++;
		next;
	}
	if ($line =~ m/Cnidaria/) {
#		print "Found Cnidaria\n";
		$flag = 0;
		$i++;
		next;
	}
	if (($flag==1) & ($line =~ m/species$/)) {
#		print "Found Chordata species\n";
		@line = split(/\*/, $line);
		$species = $line[1];
		$rdp{$species} = 1;
#		print $species."\n";
	}

	$i++;
	$line=();
	$species = ();
}
$i=0;

# count total number of RDP Chordata USA and Canada species
$num = keys %rdp;
print "There are ", $num, " unique species in RDP classifier.\n";
$num=();

# parse the nodes file to grab all species level taxids
while ($nodes[$i]) {
	my $line = $nodes[$i];
	chomp $line;

	@line = split(/\t\|\t/, $line);
	$taxid = $line[0];
	$rank = $line[2];
#	print "taxid: ", $taxid, "\trank: ", $rank."\n";

	if ($rank eq "species") {
#		print "nodes: got here\n";
#		@line = split('\t|\t', $line);
#		$taxid = $line[0];
		$taxids{$taxid} = 1;
	}
	$i++;
	$line=();
	$taxid=();
}
$i=0;

# parse names file to populate scientific names and synonyms hashes
while ($names[$i]) {
	my $line = $names[$i];
	chomp $line;

	@line = split(/\t\|/, $line);
	$taxid = $line[0];
#	print "...",$taxid,"...\n";
	$name_txt = $line[1];
#	print $name_txt."\n";
	$name_txt =~ s/\t//g;
#	print $name_txt."\n";
	$name_txt =~ s/ /_/g;
#	print $name_txt."\n";
	$name_class = $line[3];
	$name_class =~ s/\t//g;
#	print "...",$name_class."...\n";

	if (exists $taxids{$taxid}) {
#		print "found species level taxid\n";
		if ($name_class eq 'scientific name') {
#			print "found scientific name"."\n";
			$scientificNames{$taxid} = $name_txt;
			$scientificNames2{$name_txt} = $taxid;
#			print $name_txt."\n";
		}

		elsif ($name_class eq 'synonym') {
#			print "found synonym\n";
			$synonyms{$name_txt} = $taxid;
		}
	}

	$i++;
	$line=();
	@line=();
	$taxid=();
	$name_txt=();
	#$unique_name=();
	$name_class=();
}
$i=0;

#print Dumper(\%synonyms);

# loop through the GBIF North American Chordata records and see how many are in RDP Classifier
while (my ($species, $value) = each %gbif) {
	if (exists $rdp{$species}) {
		$i++;
	}
	elsif (exists $synonyms{$species}) {
		$j++;
		$taxid = $synonyms{$species};
		$name_txt = $scientificNames{$taxid};
		if (exists $rdp{$name_txt}) {
			$i++;
		}
	}
	elsif (exists $scientificNames2{$species}) {
		$k++;
	}
	else {
#		print "Can't find ", $species, " in RDP classifier or even GenBank\n";
		$l++;
	}
}

print "The RDP Classifier contains ", $i, " of the GBIF species.\n";
print "This includes ", $j, " successfully mapped synonyms from GBIF to NCBI taxonomy.\n";
print "Another ", $k, " were present in the NCBI, but not in the 12S v3 classifier.\n";
print "Finally ", $l, " were not found in NCBI at all.\n";

$i=0;
$j=0;


#!/usr/bin/perl
#March 18, 2013 by Terri Porter
#Script to get number of unique genera and species represented in processed.fasta from sample_kmer2.plx
#usage perl get_uniq_genera_species_processed_fasta.plx processed.fasta gb_org.map.uniq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb;
my $species;
my $genus;
my $uniq_genus;
my $uniq_species;
my $map;

#declare array
my @fasta;
my @map;
my @species;
my @mapline;

#declare hash
my %genus; #indexed by gb
my %species; #indexed by gb
my %uniq_genus;
my %uniq_species;

open (FASTA, "<", $ARGV[0]) || die "Error cannot open processed.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

open (MAP, "<", $ARGV[1]) || die "Error cannot open gb_org.map.uniq: $!\n";
@map = <MAP>;
close MAP;

#add genus and species names to hash indexed by gb accession
while ($map[$i]) {
	$map = $map[$i];
	chomp $map;

	@mapline = split(/\t/, $map);
	$gb = $mapline[0];
	$species = $mapline[1];
	@species = split(" ", $species);
	$genus = $species[0];

	print "MAP genus:$genus\tspecies:$species\n";

	$genus{$gb} = $genus;
	$species{$gb} = $species;
		
	$i++;
	$map=();
	@mapline=();
	$gb=();
	$species=();
	@species=();
	$genus=();
}
$i=0;

#grab genus and species for each gb in processed.fasta and add to hashes
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$gb = $line;
		$gb =~ s/^>//;
		print "FASTA gb:$gb\n";

		if (exists $genus{$gb}) {
			$genus = $genus{$gb};
			$uniq_genus{$genus} = 1;
		}

		if (exists $species{$gb}) {
			$species = $species{$gb};
			$uniq_species{$species} = 1;
		}

	}

	$i++;
	$line=();
	$gb=();
}
$i=0;

#count number of unique genera and species

$uniq_genus = scalar keys %uniq_genus;

$uniq_species = scalar keys %uniq_species;

print "unique genera: $uniq_genus\nunique species: $uniq_species\n";

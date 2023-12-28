#!/usr/bin/perl
#June 7, 2012 by Terri Porter
#Script to compare gb_org.map (from quick_filter_gb_results.plx) with readid_taxonid.parsed (from MEGAN, taxonomy_crawl.plx, and taxonomy_crawl_map.plx)
#usage perl compare_taxonomy_maps.plx readid_taxonid.parsed gb_org.map

use warnings;
use strict;

#declare var
my $i=0;
my $line;
my $asstId;
my $asstGB;
my $asstGenus;
my $refGB;
my $refSpecies;
my $refGenus;
my $match=0;
my $mismatch=0;
my $total=0;
my $original;
my $new;
my $flag=0;
my $baitMATCH=0;
my $baitMISMATCH=0;
my $baitTOTAL=0;
my $baitNotGenus=0;
my $asstSpecies; #new
my $baitNotSpecies=0;

#declare array
my @megan;
my @ref;
my @asstId;
my @line;
my @refSpecies;
my @asstGenus;
my @asstSpecies;

#declare hash
my %asst;
my %ref;
my %asstSpecies;
my %refSpecies;

open (MEGAN, "<", $ARGV[0]) || die "Error cannot open parsed MEGAN infile: $!\n";
@megan = <MEGAN>;
close MEGAN;

open (REF, "<", $ARGV[1]) || die "Error cannot open reference infile: $!\n";
@ref = <REF>;
close REF;

#add MEGAN id and genus to hashes for genera and for species
while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$asstId = $line[0];
	$asstGenus = $line[7];
	$asstSpecies = $line[8]; #genus and species here
#	@asstSpecies = split(" ",$asstSpecies);
#	$asstSpecies = $asstSpecies[1];

	if ($asstGenus ne "") {
		if (exists $asst{$asstGenus}) {
			$original = $asst{$asstGenus};
			$new = $original."|".$asstId;
			$asst{$asstGenus} = $new;
		}
		else {
			$asst{$asstGenus} = $asstId;
		}
	}
	
	print "ASST:$asstId\t$asstGenus\t$asstSpecies\n";

	if ($asstSpecies ne "") {
		if (exists $asstSpecies{$asstSpecies}) {
			$original = $asstSpecies{$asstSpecies};
			$new = $original."|".$asstId;
			$asstSpecies{$asstSpecies} = $new;
		}
		else {
			$asstSpecies{$asstSpecies} = $asstId;
		}
	}

	$i++;
	$line=();
	@line=();
	$asstId=();
	@asstId=();
	$asstGB=();
	$asstGenus=();
	$original=();
	$new=();
	$asstSpecies=();
	@asstSpecies=();
}
$i=0;

#add ref id and genus and species to hashes
while ($ref[$i]) {
	$line = $ref[$i];
	chomp $line;
	
	if ($line ne "") {
		@line = split(/\t/, $line);
		$refGB = $line[0];
		$refSpecies = $line[1];
		@refSpecies = split(" ", $refSpecies);
		$refGenus = $refSpecies[0];
		print "REF:$refGB\t$refGenus\t$refSpecies\n";
		$ref{$refGenus} = $refGB;
		$refSpecies{$refSpecies} = $refGB;
	}

	$i++;
	$line=();
	@line=();
	$refGB=();
	$refSpecies=();
	@refSpecies=();
	$refGenus=();
	$refSpecies=();
}
$i=0;

#compare ref with MEGAN asst

open (OUT, ">>", "compare_genus.txt") || die "Error cannot open compare.genus outfile: $!\n";

while (($refGenus,$refGB) = each (%ref)) {

	if (exists $asst{$refGenus}) {
		print OUT "$refGenus\t$asstGenus\tMATCH\n";
		$match++;
	}
	else {
		print OUT "$refGenus\t\tNOTFOUND\n";
		$mismatch++;
	}
	$total++;
	$flag=0;
	@asstGenus=();
	$refGenus=();
}

print "GENUS REFERENCE:\n\n$match MATCH\n$mismatch NOTFOUND\n$total TOTAL\n\n";
close OUT;
$match=0;
$flag=0;
$total=0;
$mismatch=0;

open (OUT2, ">>", "compare_species.txt") || die "Error cannot open compare.species outfile: $!\n";

while (($refSpecies,$refGB) = each (%refSpecies)) {

	if (exists $asstSpecies{$refSpecies}) {
		print OUT2 "$refSpecies\t$asstSpecies\tMATCH\n";
		$match++;
	}
	else {
		print OUT2 "$refSpecies\t\tNOTFOUND\n";
		$mismatch++;
	}
	$total++;
	$flag=0;
	@asstSpecies=();
	$refSpecies=();
}

print "SPECIES REFERENCE:\n\n$match MATCH\n$mismatch NOT FOUND\n$total TOTAL\n\n";
close OUT2;

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

#declare array
my @megan;
my @ref;
my @asstId;
my @line;
my @refSpecies;
my @asstGenus;

#declare hash
my %asst;
my %ref;

open (MEGAN, "<", $ARGV[0]) || die "Error cannot open parsed MEGAN infile: $!\n";
@megan = <MEGAN>;
close MEGAN;

open (REF, "<", $ARGV[1]) || die "Error cannot open reference infile: $!\n";
@ref = <REF>;
close REF;

#add MEGAN id and genus to hash
while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$asstId = $line[0];
	@asstId = split("_", $asstId);
	$asstGB = $asstId[0];
	$asstGenus = $line[7];

	if (exists $asst{$asstGB}) {
		$original = $asst{$asstGB};
		$new = $original."|".$asstGenus;
		$asst{$asstGB} = $new;
	}
	else {
		$asst{$asstGB} = $asstGenus;
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
}
$i=0;

#add ref id and genus to hash
while ($ref[$i]) {
	$line = $ref[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$refGB = $line[0];
	$refSpecies = $line[1];
	@refSpecies = split(" ", $refSpecies);
	$refGenus = $refSpecies[0];
	$ref{$refGB} = $refGenus;

	$i++;
	$line=();
	@line=();
	$refGB=();
	$refSpecies=();
	@refSpecies=();
	$refGenus=();
}
$i=0;

#compare ref with MEGAN asst

open (OUT, ">>", "compare.txt") || die "Error cannot open compare.txt outfile: $!\n";

while (($asstGB,$asstGenus) = each (%asst)) {

	@asstGenus = split(/\|/, $asstGenus);
	
	foreach $asstGenus (@asstGenus) {
		if (exists $ref{$asstGB}) {
			$refGenus = $ref{$asstGB};
			if ($asstGenus eq $refGenus) {
				$flag = 1;
				$baitMATCH++;
			}
			elsif ($asstGenus eq "undef") {
				$baitNotGenus++;
			}
			elsif ($asstGenus ne $refGenus) {
				$baitMISMATCH++;
			}
		}
		$baitTOTAL++;
	}

	if ($flag == 1) {
		print OUT "$asstGB\t$refGenus\t$asstGenus\tMATCH\n";
		$match++;
	}
	else {
		print OUT "$asstGB\t$refGenus\t$asstGenus\tMISMATCH\n";
		$mismatch++;
	}
	$total++;
	$flag=0;
	@asstGenus=();
	$refGenus=();
}

print "REFERENCE:\n\n$match MATCH\n\n";
print "BAIT:\n\n$baitMATCH MATCH\n$baitNotGenus NOTGENUS\n$baitMISMATCH MISMATCH\n$baitTOTAL TOTAL\n\n";

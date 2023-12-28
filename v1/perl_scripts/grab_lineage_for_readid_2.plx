#!/usr/bin/perl
#Edited Aug. 29, 2013 to start with mothurIDS (jsut the number without the Otu000 prefix)
#Aug. 27, 2013 by Terri Porter
#Script to accept list of readids from file and provide corresponding lineage for centroid
#make sure that readid_OTUcentroid.map and readid_taxonid.parsed are available
#usage perl grab_lineage_for_readid.plx mothurids.txt mothurOTU_readid.map readid_OTUcentroid.map readid_taxonid.parsed

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $readid;
my $centroid;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $lineage;
my $mothurid;

#declare array
my @map;
my @line;
my @in;
my @lineage;
my @mothur;

#declare hash
my %readid_centroid; #indexed by readid
my %lineage; #indexed by centroid
my %mothur; #indexed by mothurid (just the number wihtout the Otu000 prefix

#hash mothurids.txt
open (MOTHUR, "<", $ARGV[1]) || die "Error cannot open mothur ids: $!\n";
@mothur = <MOTHUR>;
close MOTHUR;

while ($mothur[$i]) {
	$line = $mothur[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$mothurid = $line[0];
	$readid = $line[1];

	$mothur{$mothurid} = $readid;

	$i++;
	$line=();
	@line=();
	$mothurid=();
	$readid=();
}
$i=0;

#hash readid_OTUcentroid.map
open (MAP, "<", $ARGV[2]) || die "Error cannot open map: $!\n";
@map = <MAP>;
close MAP;

while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$readid = $line[0];
	$centroid = $line[1];

	$readid_centroid{$readid} = $centroid;

	$i++;
	$line=();
	@line=();
	$readid=();
	$centroid=();
}
$i=0;

#hash readid_taxonid.parsed
open (LINEAGE, "<", $ARGV[3]) || die "Error cannot open lineage file: $!\n";
@lineage = <LINEAGE>;
close LINEAGE;

while ($lineage[$i]) {
	$line = $lineage[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$centroid = $line[0];
	$kingdom = $line[2];
	$phylum = $line[3];
	$class = $line[4];
	$order = $line[5];
	$family = $line[6];
	$genus = $line[7];

	$lineage{$centroid} = $kingdom."|".$phylum."|".$class."|".$order."|".$family."|".$genus;

	$i++;
	$line=();
	@line=();
	$centroid=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
}
$i=0;

#print it out!
open (IN, "<", $ARGV[0]) || die "Error cannot open mothurids.txt:$!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "indicator_species_lineages.txt") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$mothurid = $in[$i];
	chomp $mothurid;

	if (exists $mothur{$mothurid}) {
		$readid = $mothur{$mothurid};
		
		if (exists $readid_centroid{$readid}) {
			$centroid = $readid_centroid{$readid};
	
			if (exists $lineage{$centroid}) {
				$lineage = $lineage{$centroid};
				@lineage = split(/\|/,$lineage);
				$kingdom = $lineage[0];
				$phylum = $lineage[1];
				$class = $lineage[2];
				$order = $lineage[3];
				$family = $lineage[4];
				$genus = $lineage[5];

				print OUT "$mothurid\t$readid\t$centroid\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\n";

			}
			else {
				print "Warning: can't find lineage for $centroid\n";
			}
		}
		else {
			print "Warning: can't find centroid for $readid\n";
		}
	}
	else {
		print "Warning: can't find readid for $mothurid\n";
	}
	
	$i++;
	$mothurid=();
	$readid=();
	$centroid=();
	$lineage=();
	@lineage=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
}
$i=0;

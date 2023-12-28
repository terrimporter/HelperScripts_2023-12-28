#!/usr/bin/perl
#Aug. 27, 2013 by Terri Porter
#Script to accept list of readids from file and provide corresponding lineage for centroid
#make sure that readid_OTUcentroid.map and readid_taxonid.parsed are available
#usage perl grab_lineage_for_readid.plx readids.txt readid_OTUcentroid.map readid_taxonid.parsed

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

#declare array
my @map;
my @line;
my @in;
my @lineage;

#declare hash
my %readid_centroid; #indexed by readid
my %lineage; #indexed by centroid

#hash readid_OTUcentroid.map
open (MAP, "<", $ARGV[1]) || die "Error cannot open map: $!\n";
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
open (LINEAGE, "<", $ARGV[2]) || die "Error cannot open lineage file: $!\n";
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
open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$readid = $in[$i];
	chomp $readid;

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

			print "$readid\t$centroid\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\n";

		}
		else {
			print "Warning: can't find lineage for $centroid\n";
		}
	}
	else {
		print "Warning: can't find centroid for $readid\n";
	}
	
	$i++;
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

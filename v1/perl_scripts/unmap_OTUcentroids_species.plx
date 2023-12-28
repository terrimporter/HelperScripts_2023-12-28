#!/usr/bin/perl
#Sept.13/13 add species rank
#Terri Porter, May 30, 2013
#Script to unmap OTUcentroids after MEGAN processing, taxonoy_crawl.plx, and taxonomy_cral_map.plx
#Use to grab samples related to each OTUcentroid classified so that it can be sorted in excel
#usage perl unmap_OTUcentroids.plx readid_taxonid.parsed readid_OTUcentroid.map readid_sample.map

use strict;
use warnings;

#var
my $i=0;
my $line;
my $readid;
my $OTUcentroid;
my $sample;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $lineage;
my $taxonid;
my $detail;

#array
my @in;
my @map1;
my @map2;
my @line;
my @in2;
my @detail;

#hash
my %OTUcentroid; #indexed by readid
my %sample; #indexed by readid, then cleared, then later indexed by sample
my %detail; #indexed by readid;

open (MAP1, "<", $ARGV[1]) || die "Error can't open first mapping file: $!\n";
@map1 = <MAP1>;
close MAP1;

#hash map1
while ($map1[$i]) {
	$line = $map1[$i];
	chomp $line;

	if (length ($line) > 0 ) {
		@line = split(/\t/, $line);
		$readid = $line[0];
		$OTUcentroid = $line[1];
		$OTUcentroid{$readid} = $OTUcentroid;
	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$OTUcentroid=();
}
$i=0;

open (MAP2, "<", $ARGV[2]) || die "Error can't open second mapping file: $!\n";
@map2 = <MAP2>;
close MAP2;

#hash map2
while ($map2[$i]) {
	$line = $map2[$i];
	chomp $line;

	if (length ($line) > 0 ) {
		@line = split(/\t/, $line);
		$readid = $line[0];
		$sample = $line[1];
		$sample{$readid} = $sample;

	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$sample=();
}
$i=0;

open (IN, "<", $ARGV[0]) || die "Error can't open first infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "unmapped_readids.txt") || die "Error can't open outfile: $!\n";

open (OUT2, ">>", "unmapped_OTUcentroids.txt") || die "Error can't open second outfile: $!\n";

print OUT "Sample\tOTUcentroid\tReadid\tTaxonid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
print OUT2 "Sample\tOTUcentroid\tTaxonid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";

#map readid to OTUcentroids, retaining OTUcentroid lineage for each
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if (length ($line) > 0 ) {
		@line = split(/\t/, $line);
		$OTUcentroid = $line[0];
		$OTUcentroid =~ s/' '//g; #remove any hidden spaces
		$taxonid = $line[1];
		$kingdom = $line[2];
		$phylum = $line[3];
		$class = $line[4];
		$order = $line[5];
		$family = $line[6];
		$genus = $line[7];
		$species = $line[8];
		$lineage = $kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species;

		while ( my ($key,$val) = each (%OTUcentroid) ) {
			$val =~ s/' '//g; #remove any hidden spaces
			if ( $val eq $OTUcentroid) {
				if ( exists $sample{$val}) {
					$sample = $sample{$val};
					print OUT $sample."\t".$val."\t".$key."\t".$taxonid."\t".$lineage."\n";
				}
				else {
					print "Can't find sample for read $val.\n";
				}
			}
		}
	}
	$i++;
	$line=();
	@line=();
	$OTUcentroid=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
	$lineage=();
}
$i=0;
close OUT;

%sample=();

open (IN2, "<", "unmapped_readids.txt") || die "Error can't open unmapped_readids.txt: $!\n";
@in2 = <IN2>;
close IN2;

#hash file to get list of unique samples, and details
while ($in2[$i]) {
	$line = $in2[$i];
	chomp $line;

	if ($line =~ /^Sample/) {
		$i++;
		next;
	}

	else {
		@line = split(/\t/, $line);
		$sample = $line[0];
		$OTUcentroid = $line[1];
		$readid = $line[2];
		$taxonid = $line[3];
		$kingdom = $line[4];
		$phylum = $line[5];
		$class = $line[6];
		$order = $line[7];
		$family = $line[8];
		$genus = $line[9];
		$species = $line[10];
		$lineage = $taxonid."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species;
		$detail = $sample."\|".$OTUcentroid;
#		$sample{$sample} = 1;  #just to get a list of unique samples ###or provide a list of samples to parse for!
		$detail{$detail} = $lineage;
	}
	$i++;
	$line=();
	@line=();
}
$i=0;

##### custom list of samples for Mehrdad #####

%sample=();
$sample{'PAD'} = 1;
$sample{'WC'} = 1;

#####

#print out OTUcentroids for each sample
while (my ($key,$val) = each(%sample) ) {
		
	while (($detail,$lineage) = each(%detail) ) {
		@detail = split(/\|/,$detail);
		$sample = $detail[0];
		$OTUcentroid = $detail[1];
		
#		print "$key\t$sample\n";

		if ($sample  =~ /$key/) {
			print OUT2 $sample."\t".$OTUcentroid."\t".$lineage."\n";
		}
	}
}
close OUT2;

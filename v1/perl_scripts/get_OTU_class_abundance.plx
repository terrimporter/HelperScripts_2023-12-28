#!/usr/bin/perl
#April 17, 2013 by Terri Porter
#Script to use MEGAN csv readid-taxonid output at the 'class' rank and 
#Usearch otus_minsize2.fa.newlineremoved (from run_usearch_OTU_clustering.plx and remove_newline.plx) 
#to get OTUreadid\ttaxonid\tOTUreadidAbundance
#Manually manipulate this file in excel to generate data for Nagissa's100% stacked bar graphs
#USAGE perl get_OTU_class_abundance.plx MEGAN-ex.txt otus_minsize2.fa.newlineremoved

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $readid;
my $taxonid;
my $size;

#declare array
my @megan;
my @usearch;

#declare hash
my %taxonid; #indexed by readid
my %size; #indexed by readid

open (MEGAN, "<", $ARGV[0]) || die "Error cannot open MEGAN file: $!\n";
@megan = <MEGAN>;
close MEGAN;

open (USEARCH, "<", $ARGV[1]) || die "Error cannot open USEARCH file: $!\n";
@usearch = <USEARCH>;
close USEARCH;

#hash OTUreadids and taxonids
while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;

	if ($line !~ /^\s+/) {
		$line =~ /^(\w+)/;
		$readid = $1;
		$line =~ /(\d+)$/;
		$taxonid = $1;
		$taxonid{$readid} = $taxonid;
	}
	$i++;
	$line=();
	$readid=();
	$taxonid=();
}
$i=0;

#hash OTUreadids and cluster size
while ($usearch[$i]) {
	$line = $usearch[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w{14})/;
		$readid = $1;
		$line =~ /size=(\d+)/;
		$size = $1;
#		print "size:$size\n";#test
		$size{$readid} = $size;
	}
	$i++;
	$line=();
	$readid=();
	$size=();
}
$i=0;

#put the data together into a table for excel
open (OUT, ">>", "OTU_abundance.txt") || die "Error cannot open outfile: $!\n";

while ( ($readid,$taxonid) = each (%taxonid) ) {
	if (exists $size{$readid}) {
		$size = $size{$readid};
		print OUT "$taxonid\t$readid\t$size\n";
	}
	else {
		print "Error cannot find OTU size for $readid\n";
	}
	$size=();
}
close OUT;

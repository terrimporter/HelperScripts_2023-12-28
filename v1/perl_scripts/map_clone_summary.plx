#!/usr/bin/perl
#Oct.6, 2010 by Terri Porter
#Script to swap out accession id for genus and species.  Do this prior to turning this into a table for heatmap in R.
#usage $perl map_clone_summary.plx clone_summary.txt original_genbank.map

use strict;
use warnings;

#declare var
my $line;
my $accession_full;
my $accession_short;
my $genus;
my $species;
my $i=0;
my $j=0;
my $line2;
my $accession2;
my $genus2;
my $species2;
my $pp2;

#declare array
my @clone_summary;
my @original_genbank;
my @line;
my @accession_full;
my @line2;

open (IN1,"<",$ARGV[0]) || die ("Error:$!\n");
@clone_summary = <IN1>;
close IN1;

open (IN2,"<",$ARGV[1]) || die ("Error:$!\n");
@original_genbank = <IN2>;
close IN2;

open (OUT,">>","clone_summary.mapped") || die ("Error:$!\n");

while ($original_genbank[$i]) {
	$line = $original_genbank[$i];
	chomp $line;
	@line = split ( /\t/,$line);
	$accession_full = $line[0];
	if ($accession_full =~ /\./) {
		@accession_full = split(/\./,$accession_full);
		$accession_short = $accession_full[0];
	}
	else {
		$accession_short = $accession_full;
	}
	$genus = $line[1];
	$species = $line[2];
	
	while ($clone_summary[$j]) {
		$line2 = $clone_summary[$j];
		chomp $line2;
		if ($line2 =~ /$accession_short/) {
			@line2 = split (/\t/,$line2);	
			$accession2 = $line2[0];
			$genus2 = $line2[1];
			$species2 = $line2[2];
			$pp2 = $line2[3];
			print OUT "$accession_short\t$genus\t$species\t$genus2\t$species2\t$pp2\n";
			$j++;
		}
		else {
			$j++;
			next;
		}
	}
	$j=0;
	$i++;
}

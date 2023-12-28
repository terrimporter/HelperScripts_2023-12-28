#!/usr/bin/perl
#June 28, 2013 by Terri Porter
#Script to parse unmapped_OTUcentroids.txt from unmap_OTUcentroids.plx to create table for stacked bar charts in excel
#usage perl parse_unampped_OTUcentroids.plx unmapped_OTUcentroids.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $sample;
my $family;
my $original;
my $pad;
my $wc;
my $string;
my $taxonid;

#declare array
my @in;
my @line;
my @original;
my @string;

#declare hash
my %family; #indexed by family name

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>; 
close IN;

#parse and hash the infile
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($i > 0) {
		@line = split(/\t/, $line);
		$sample = $line[0];

		if ($sample =~ /^(F|R)_/) { #only for pooled 16S v3 and v6 sets
			$sample =~ s/^(F|R)_//g;
		}

		if ($sample =~ /^PAD/) {
			$sample = "PAD";
		}
		elsif ($sample =~ /^WC/) {
			$sample = "WC";
		}

		$taxonid = $line[2];

		if ($taxonid == "-1") {
			$family = "No BLAST hits";
		}
		elsif ($taxonid == "-2") {
			$family = "Not assigned";
		}
		elsif ($taxonid == "311242") { #ITS4F
			$family = "311242_1302847";
		}
		elsif ($taxonid == "70108") {
			$family = "70108_1295096";
		}
		elsif ($taxonid == "387323") {
			$family = "387323_1077376";
		}
		elsif ($taxonid == "470448") { #ITS1R
			$family = "470447_1307819";
		}
		elsif ($taxonid == "43432") {
			$family = "43432_1250544";
		}
		elsif ($taxonid == "480077") {
			$family = "480077_1075807";
		}
		elsif ($taxonid == "4889") {
			$family = "4889_1213189";
		}
		elsif ($taxonid == "475271") {
			$family = "475271_176275";
		}
		elsif ($taxonid == "421611") { #16S v6
			$family = "421611_74719";
		}
		else {
			$family = $line[8]; #order [8]; family [7]; genus [6]
			if (length($family) == 0) {
				print "problem parsing taxonid: $taxonid\n";
			}
			elsif (length($family) eq '') {
				print "problem parsing taxonid: $taxonid\n";
			}
		}

		if (exists $family{$family}) {
			$original = $family{$family};
			@original = split(/\|/, $original);
			$pad = $original[0];
			$wc = $original[1];

			if ($sample eq "PAD") {
				$pad++;
			}
			elsif ($sample eq "WC") {
				$wc++;
			}
			$family{$family} = $pad."|".$wc;
		}
		else {
			if ($sample eq "PAD") {
				$family{$family} = "1"."|"."0";
			}
			elsif ($sample eq "WC") {
				$family{$family} = "0"."|"."1";
			}
		}
	}

	$i++;
	$line=();
	@line=();
	$sample=();
	$family=();
	$original=();
	@original=();
	$pad=();
	$wc=();
}
$i=0;

#print out the sorted hash
foreach $family (sort keys %family) {
	$string = $family{$family};
	@string = split(/\|/, $string);
	$string = join("\t", @string);
	print "$family\t$string\n";
}


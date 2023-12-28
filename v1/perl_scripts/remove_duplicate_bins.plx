#!/usr/bin/perl
#May 14, 2012 by Terri Porter
#Script to remove BINs from bin.query that are already in bin.html dir
#in html dir/, ls > html.list, sed 's/\.html//g' html.list > html.bin
#usage perl remove_duplicate_bins.plx BIN.query html.bin

use strict;
use warnings;

#declare vary
my $i=0;
my $html;
my $line;
my $j=0;
my $BIN;

#declare array
my @query;
my @duplicates;
my @line;

#declare hash
my %html;

open (IN, "<", $ARGV[0]) || die "Error cannot open first infile: $!\n";
@query = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Error cannot open second infile: $!\n";
@duplicates = <IN2>;
close IN2;

while ($duplicates[$i]) {
	$html = $duplicates[$i];
	chomp $html;
	$html{$html} = 1;
	$i++;
}
$i=0;

open (OUT, ">>", "NEW_BIN.query") || die "Error cannot open outfile: $!\n";

while ($query[$i]) {
	$line = $query[$i];
	chomp $line;

	@line = split(" ", $line);

	while ($line[$j]) {
		$BIN = $line[$j];

		if (exists $html{$BIN}) {
			$j++;
			next;
		}
		else {
			print OUT "$BIN ";
		}
		$j++;
		$BIN=();
	}
	$j=0;

	$i++;
	$line=();
	@line=();
}
$i=0;
close OUT;

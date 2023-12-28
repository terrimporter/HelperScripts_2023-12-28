#!/usr/bin/perl
#March 21, 2013 by Terri Porter
#Script to grab gb accessions from gb_seq.map.uniq from quick_filter_gb_results.plx and use these to grab associated taxids from gb_taxid.map.uniq
#usage perl grab_taxids_from_gb.plx gb_seq.map.uniq gb_taxid.map.uniq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb;
my $tax;
my $seq;

#declare array
my @seq;
my @tax;
my @line;
my @gb;

#declare hash
my %gb; #indexed by gb
my %tax; #indexed by gb

open (SEQ, "<", $ARGV[0]) || die "Error cannot open gb_seq.map.uniq: $!\n";
@seq = <SEQ>;
close SEQ;

open (TAX, "<", $ARGV[1]) || die "Error cannot open gb_taxid.map.uniq: $!\n";
@tax = <TAX>;
close TAX;

#hash gb_seq file
while ($seq[$i]) {
	$line = $seq[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/, $line);
		$gb = $line[0];
		$seq = $line[1];
		$gb{$gb} = $seq;
	}
	$i++;
	$line=();
	@line=();
	$gb=();
	$seq=();
}
$i=0;

#hash gb_taxid file
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/, $line);
		$gb = $line[0];
		$tax = $line[1];
		$tax{$gb} = $tax;
	}
	$i++;
	$line=();
	@line=();
	$gb=();
	$tax=();
}
$i=0;

#grab just the taxids I need and print out a list for taxonomy_crawl_species.plx
open (OUT, ">>", "filtered_taxids.txt") || die "Error cannot open outfile: $!\n";

while (my ($gb, $seq) = each %gb) {

	if (exists $tax{$gb}) {
		$tax = $tax{$gb};
		print OUT "$tax\n";
	}
	else {
		print "Taxid for $gb not found\n";
	}
}
close OUT;

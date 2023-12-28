#!/usr/bin/perl
# July 4, 2020 based on grab_taxids_from_gb_2.plx
#April 10, 2018 add a step to check sequence for non-nucleotide characters (anything other than A,C,G,T)
#March 21, 2013 by Terri Porter
#Script to grab gb accessions from gb_seq.map.uniq from quick_filter_gb_results.plx and use these to grab associated taxids from gb_taxid.map.uniq
#usage perl filter_taxids_and_seqs.plx gb_seq.map.uniq gb_taxid.map.uniq

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
my %gb_poorqual; #indexed by gb
my %filtered_taxids; # key=taxid, value = 1

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

		# Check for non-nucleotide characters here (anything other than an A, C, T, or G)
		if ($seq =~ /[^ACGTacgt]/ ) {
			$gb_poorqual{$gb} = 1;
			$i++;
			next;
		}
		else {
			$gb{$gb} = $seq;
		}
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
		if (exists $gb_poorqual{$gb}) {
			$i++;
			next;
		}
		else {
			$tax{$gb} = $tax;
		}
	}
	$i++;
	$line=();
	@line=();
	$gb=();
	$tax=();
}
$i=0;

#grab just the taxids I need and print out a list for taxonomy_crawl_species.plx
while (my ($gb, $seq) = each %gb) {

	if (exists $tax{$gb}) {
		$tax = $tax{$gb};
		$filtered_taxids{$tax} = 1;
	}
}

open (OUT, ">>", "filtered_taxids.txt") || die "Error cannot open outfile: $!\n";

while (my ($taxid, $value) = each %filtered_taxids) {
	print OUT "$taxid\n";
}
close OUT;

#print out a file of poor quality GBs to avoid in next processing steps
open (OUT2, ">>", "gb_poorqual.txt") || die "Error cannot open outfile2: $!\n";

while (my ($gb, $null) = each %gb_poorqual) {

	print OUT2 "$gb\n";

}
close OUT2;

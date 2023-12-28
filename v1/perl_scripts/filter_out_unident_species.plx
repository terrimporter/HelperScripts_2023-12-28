#!/usr/bin/perl
# Teresita M. Porter, Feb. 8, 2021
# Script to exclude records unidentified at species rank
# USAGE perl filter_out_unident_species.plx unite_outgroup.fasta unite_outgroup.txt

use strict;
use warnings;
use Data::Dumper;

# declare var
my $i=0;
my $line;
my $sh;
my $lineage;
my $speciesField;
my $fas_out = "unite_outgroup2.fasta";
my $tax_out = "unite_outgroup2.txt";
my $j;
my $seq;

# declare array
my @fas;
my @tax;
my @line;
my @lineage;

# declare hashes
my %tax; # key = sh, val = lineage

open (TAX, "<", $ARGV[1]) || die "Error cannot open taxonomy file: $!\n";
@tax = <TAX>;
close TAX;

# hash taxonomy
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$sh = $line[0];
	$lineage = $line[1];

	@lineage = split(/;/, $lineage);
	$speciesField = pop @lineage;

	if ($speciesField =~ /unidentified/) {
		$i++;
		next;
	}
	else {
		$tax{$sh} = $lineage;
	}

	$i++;

}
$i=0;


open (FAS, "<", $ARGV[0]) || die "Error cannot open fasta file: $!\n";
@fas = <FAS>;
close FAS;

# for each sh in hash, print out new fasta and taxonomy files
open (FAS2, ">>", $fas_out) || die "Error cannot open fasta outfile: $!\n";
open (TAX2, ">>", $tax_out) || die "Error cannot open taxonomy outfile: $!\n";

while($fas[$i]) {
	$line = $fas[$i];
	chomp $line;

	if ($line =~ /^>/) { # header
		@line = split(/\t/,$line);
		$sh = $line[0];
		$sh =~ s/^>//;

		if (exists $tax{$sh}) {
			print FAS2 $line."\n";
			$j = $i+1;
			$seq = $fas[$j];
			chomp $seq;
			print FAS2 $seq."\n";
			$lineage = $tax{$sh};
			print TAX2 "$sh\t$lineage\n";
		}
		else {
			$i+=2;
			next;
		}
	}
	$i++;
}
$i=0;
close FAS2;
close TAX2;

#!/usr/bin/perl
# Teresita M. Porter, Jan. 23/21
# Script to remove singletons and doubletons from cat.denoised.nonchimeras file
# USAGE perl filter_denoised_nonchimeras.plx ESV.table.parsed cat.denoised.nonchimeras

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $zotu;
my $outfile="cat.denoised.nonchimeras.parsed";
my $flag = 0;

# declare arrays
my @table;
my @line;
my @fasta;

# declare hashes
my %zotus; # key = zotu, value = 1

open (TABLE, "<", $ARGV[0]) || die "Error cann't open ESv.table.parsed: $!\n";
@table = <TABLE>;
close TABLE;

# hash Zotus to keep (singletons and doubletons already removed)
while ($table[$i]) {
	$line = $table[$i];
	chomp $line;

	if ($i == 0) {
		$i++;
		next;
	}
	else {
		@line = split(/\t/, $line);
		$zotu = shift @line;
		$zotus{$zotu} = 1;
	}

	$i++;

}
$i=0;

open (FASTA, "<", $ARGV[1]) || die "Error canot open cat.denoised.nonchimeras: $!\n";
@fasta = <FASTA>;
close FASTA;

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) { # header
		$zotu = $line;
		$zotu =~ s/^>//g;

		if (exists $zotus{$zotu}) {
			print OUT $line."\n";
			$flag = 1;
		}
		else {
			$flag = 0;
		}
	}
	else {
		if ($flag == 1) {
			print OUT $line."\n";
		}
	}

	$i++;

}
$i=0;

close OUT;


#!/usr/bin/perl
# Teresita M. Porter, July 27, 2021
# Script to grab gbs from downsampled FASTA file from vsearch and filter seq and taxid mapping files accordingly
# USAGE perl fasta_to_maps.plx map.fasta.clustered80 gb_seq_nonvert.uniq gb_taxid_nonvert.uniq

use strict;
use warnings;

# declare vars
my $i=0;
my $line;
my $gb;
my $outfile1 = "gb_seq_nonvert.uniq2";
my $outfile2 = "gb_taxid_nonvert.uniq2";

# declare arrays
my @fasta;
my @gbseq;
my @line;
my @gbtaxid;

# declare hashes
my %gbs; # key = gb, value = 1

# read in FASTA of downsampled outgroup seqs from vsearch
open (IN, "<", $ARGV[0]) || die "Error cannot open FASTA file: $!\n";
@fasta = <IN>;
close IN;

# grab gb's from FASTA header
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) { # header
		$gb = $line;
		$gb =~ s/^>//g;
		$gbs{$gb} = 1;

	}
	$i++;

}
$i=0;

# read in original gb-seq mapping file
open (IN2, "<", $ARGV[1]) || die "Error cannot open gb-seq mapping file: $!\n";
@gbseq = <IN2>;
close IN2;

open (OUT, ">>", $outfile1) || die "Error cannot open new gb-seq mapping file: $!\n";

# filter gb-seq mapping file
while ($gbseq[$i]) {
	$line = $gbseq[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$gb = $line[0];

	if (exists $gbs{$gb}) {
		print OUT $line."\n";
	}
	$i++;

}
$i=0;

# read in original gb-taxid mapping file
open (IN3, "<", $ARGV[2]) || die "Error cannot open original gb-taxid mapping file: $!\n";
@gbtaxid = <IN3>;
close IN3;

open (OUT2, ">>", $outfile2) || die "Error cannot open new gb-taxid mapping file: $!\n";

# filter gb-taxid mapping file
while ($gbtaxid[$i]) {
	$line = $gbtaxid[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$gb = $line[0];

	if ($gbs{$gb}) {
		print OUT2 $line."\n";
	}

	$i++;

}
$i=0;

close OUT;
close OUT2;

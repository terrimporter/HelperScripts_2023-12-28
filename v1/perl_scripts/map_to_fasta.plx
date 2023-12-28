#!/usr/bin/perl
# Teresita M. Porter, July 27, 2021
# Script to turn two mapping files into a FASTA file (to reduce redundancy using vsearch clustering)
# USAGE perl map_to_fasta.plx gb_seq_nonvert.map

use strict;
use warnings;

# declare vars
my $i=0;
my $line;
my $gb;
my $seq;
my $outfile = "map.fasta";

# declare arrays
my @gb_seq;
my @line;

# declare hashes
my %taxids; # key - taxid, value - gb

open (IN, "<", $ARGV[0]) || die "Error cannot open mapping file: $!\n";
@gb_seq = <IN>;
close IN;

open (OUT, ">", $outfile) || die "Error cannot open outfile: $!\n";

while ($gb_seq[$i]) {
	$line = $gb_seq[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$gb = $line[0];
	$seq = $line[1];

	print OUT ">".$gb."\n".$seq."\n";
	$i++;

}
$i=0;

close OUT;

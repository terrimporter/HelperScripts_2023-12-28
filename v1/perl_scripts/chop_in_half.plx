#!/usr/bin/perl
# Teresita M. Porter, Feb. 6, 2020
# Script to chop sequences in half to reduce their length for testing
# expect a STRICT fasta file
# USAGE perl chop_in_half.plx mutated_0.fasta

use strict;
use warnings;
use Data::Dumper;
use Math::Round;

# vars
my $i=0;
my $line;
my $header;
my $j;
my $seq;
my $length;
my $firstHalfLength;
my $secondHalfLength;
my $firstHalf;
my $secondHalf;
my $outfile1 = "firstHalf.fasta";
my $outfile2 = "secondHalf.fasta";

# arrays
my @in;
my @seq;

# hashes

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT1, ">>", $outfile1) || die "Cannot open outfile1: $!\n";
open (OUT2, ">>", $outfile2) || die "Cannot open outfile2: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$j = $i+1;
		$seq = $in[$j];
		chomp $seq;

		@seq = split(//, $seq);
		$length = scalar(@seq);

		# rounds to nearest integer
		$firstHalfLength = round($length/2);
		$secondHalfLength = $length - $firstHalfLength;

		$firstHalf = substr($seq, 0, $firstHalfLength);
		$secondHalf = substr($seq, $firstHalfLength, $secondHalfLength);

#		print "seq $seq\nfirst $firstHalf\tsecond $secondHalf\n";

		print OUT1 "$header\n$firstHalf\n";
		print OUT2 "$header\n$secondHalf\n";

	}
	$i++;
}
$i=0;

close OUT1;
close OUT2;


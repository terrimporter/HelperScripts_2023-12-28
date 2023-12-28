#!/usr/bin/perl
#March 15, 2013 by Terri Porter
#Script to turn kmer.freq (from sample_kmer2.plx) to kmer_uniq.fasta.
#Use kmer_uniq.fasta instead of kmer.fasta for BLAST+MEGAN against customdb from processed.fasta
#useage perl freq_to_fasta.plx kmer.freq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $bait;
my $j;

#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "kmer_uniq.fasta") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);

	$bait = $line[0];
	$j = $i+1;
	print OUT ">$j\n$bait\n";
	
	$i++;
	$line=();
	@line=();
	$bait=();
	$j=0;
}
$i=0;

close OUT;

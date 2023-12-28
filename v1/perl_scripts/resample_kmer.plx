#!/usr/bin/perl
#March 30, 2012 by Terri Porter
#Script to randomly resample kmers from kmer.fasta without replacement
#usage perl resample_kmer.plx kmer.fasta

use strict;
use warnings;
use List::Util 'shuffle';

#declare var
my $i=0;
my $line;
my $kmerID;
my $j;
my $line2;
my $kmerSeq;
my $kmerIDSeq;
my $string;
my $newstring;
my $count;
my $proportion = 0.25;#####reset this value here, 0.25, 0.50, 0.75, 1.0
my $newcount;
my $random_index;

#declare array
my @in;
my @lines;
my @keys;
my @kmerIDSeq;

#declare hash
my %kmer;
my %lines;

open (IN,"<",$ARGV[0]) || die "Error cannot read from kmer.fasta: $!\n";
@in = <IN>;
close IN;

#add kmers to hash for easy processing
while($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$kmerID = $line;

		$j=$i+1;
		$line2 = $in[$j];
		chomp $line2;
		$kmerSeq = $line2;

		$kmerIDSeq = $kmerID."|".$kmerSeq;	
		$kmer{$kmerIDSeq} = 1;
	}
	$i++;
	$line=();
	$kmerID=();
	$j=();
	$line2=();
	$kmerSeq=();
	$kmerIDSeq=();
	$string=();
	$newstring=();
}
$i=0;

#get total number of keys in hash
$count = keys(%kmer);
$newcount = int($proportion*$count);

#randomly sample keys WITHOUT REPLACEMENT up to required proportion (25%, 50%, 75%, all/100%)
open (OUT,">>","kmer_rand.fasta") || die "Error cannot write to outfile: $!\n";

while(scalar(keys(%lines)) < $newcount) {#generate a bunch of unique random index values
	$lines{int(rand($count))} = "1";
}
@lines = keys(%lines);
#print "@lines\n";#test
#my $scalar = scalar(@lines);
#print "num sampled: $scalar\n";

@keys = keys(%kmer); #put hash keys into array so they can be sampled using the random index values from above
#print "@keys\n";#test

while($lines[$i]) {
	$random_index = $lines[$i];
	$kmerIDSeq = $keys[$random_index];
	@kmerIDSeq = split(/\|/,$kmerIDSeq);
	$kmerID = $kmerIDSeq[0];
	$kmerSeq = $kmerIDSeq[1];
	print OUT ">$kmerID\n$kmerSeq\n";

	$i++;
	$random_index=();
	$kmerIDSeq=();
	@kmerIDSeq=();
	$kmerID=();
	$kmerSeq=();
}
$i=0;

	

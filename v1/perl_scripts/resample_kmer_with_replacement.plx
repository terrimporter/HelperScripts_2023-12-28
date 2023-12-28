#!/usr/bin/perl
#April 9, 2012 edited to sample WITH REPLACEMENT and grab a nubmer of kmers instead of a proportion
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
#my $proportion = 0.25;#####reset this value here, 0.25, 0.50, 0.75, 1.0
my $number = 10000;#####reset this value log10 (1K to 1million)
my $newcount;
my $random_index;
my $range;
my $scalar=0;

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
#$newcount = int($proportion*$count);

#randomly sample keys up to required proportion (25%, 50%, 75%, all/100%) NOT!
#randomly sample keys up to required NUMBER
open (OUT,">>","kmer_rand.fasta") || die "Error cannot write to outfile: $!\n";

@keys = keys(%kmer); #put hash keys into array so they can be sampled using the random index values from above
#print "@keys\n";#test

#WITH REPLACEMENT
#while(scalar(keys(%lines)) < $newcount) {
#generate a bunch of unique random index values
#	$lines{int(rand($count))} = "1";
#}
#@lines = keys(%lines);

#WITHOUT REPLACEMENT
$range = $count+1;
while($scalar < $number) {
	$random_index = int(rand($range));
	push (@lines, $random_index);
	$scalar = scalar(@lines);
	$kmerIDSeq = $keys[$random_index];
	@kmerIDSeq = split(/\|/,$kmerIDSeq);
	$kmerID = $kmerIDSeq[0];
	$kmerSeq = $kmerIDSeq[1];
	print OUT ">$i|$kmerID\n$kmerSeq\n";

	$random_index=();
	$kmerIDSeq=();
	@kmerIDSeq=();
	$kmerID=();
	$kmerSeq=();
}

#print "@lines\n";#test
#my $scalar = scalar(@lines);
#print "num sampled: $scalar\n";

@keys = keys(%kmer); #put hash keys into array so they can be sampled using the random index values from above
#print "@keys\n";#test

#while($lines[$i]) {
#	$random_index = $lines[$i];
#	print $random_index."\n";#test
#	$kmerIDSeq = $keys[$random_index];
#	@kmerIDSeq = split(/\|/,$kmerIDSeq);
#	$kmerID = $kmerIDSeq[0];
#	$kmerSeq = $kmerIDSeq[1];
#	print OUT ">$i|$kmerID\n$kmerSeq\n";

#	$i++;
#	$random_index=();
#	$kmerIDSeq=();
#	@kmerIDSeq=();
#	$kmerID=();
#	$kmerSeq=();
#}
#$i=0;

#!/usr/bin/perl
#March 28, 2012 by Terri Porter
#Script to assign kmers to all possible reference seqs, id=1.0.  Output non-assigned kmers into separate file to deal with later.
#usage perl assign_kmer.plx processed.fasta kmer.fasta

use strict;
use warnings;

#declare var
my $ref;
my $gb;
my $i=0;
my $j=0;
my $line;
my $seq;
my $kmer;
my $kmerID;
my $kmerSeq;
my $string;
my $newstring;
my $numKmer;
my $count=0;
my $length;
my $maxKmer;
my $stepSize=100;#####edit this
my $windowSize=100;#####edit this

#declare array
my @ref;
my @kmer;
my @kmerID;
my @seq;
my @string;

#declare hash
my %ref;
my %kmer;
my %match;
my %present;
my %kmer2;

open (REF,"<",$ARGV[0]) || die "Error cannot read reference.fasta: $!\n";
@ref = <REF>;
close REF;

open (KMER,"<",$ARGV[1]) || die "Error cannot read kmer.fasta: $!\n";
@kmer = <KMER>;
close KMER;

#assign each ref seq to hash
while($ref[$i]) {
	$ref = $ref[$i];
	chomp $ref;
	if ($ref =~ /^>/) {
		$ref =~ s/^>//;
		$gb = $ref;
		$j=$i+1;
		$line = $ref[$j];
		chomp $line;
		$seq = $line;
		$ref{$gb} = $seq;
		$i+=2;
	}
	else {
		$i++;
	}
}
$i=0;
$j=();
print "Reference sequences parsed\n";

#assign all kmers to hash without any duplicates
while($kmer[$i]) {
	$kmer = $kmer[$i];
	chomp $kmer;
	if ($kmer =~ /^>/) {
		$kmer =~ s/^>//;
		$kmerID = $kmer;
		$j=$i+1;
		$line = $kmer[$j];
		chomp $line;
		$kmerSeq = $line;
		if (exists($kmer{$kmerSeq})) {
			$i+=2;
			next;
		}
		else {
			$kmer{$kmerSeq} = $kmerID;
			$i+=2;
		}
	}
	else {
		$i++;
	}
}
$i=0;
$j=0;
print "Kmers parsed\n";

#create kmer mapping hash to speed up kmer matching
#don't bother to check for ambituities or min seq length here because this was done before
while(($gb,$seq) = each(%ref)) {
	@seq = split("",$seq);
	$length = scalar(@seq);
	$maxKmer = $length/$stepSize;
	$maxKmer = int($maxKmer);
	while($i<$maxKmer){
		$kmer = substr($seq,$j,$windowSize);
		if (exists($kmer2{$kmer})) {
			$string = $kmer2{$kmer}; #$gb
			$newstring = $string.",".$gb;
			$kmer2{$kmer} = $newstring;
		}
		else {
			$kmer2{$kmer} = $gb;
		}
		$i++;
		$j+=$stepSize;
		$kmer=();
	}
	$i=0;
	$j=0;
	@seq=();
	$length=();
	$maxKmer=();
}

#while(($kmer,$string) = each(%kmer2)){
#	print $kmer."\t".$string."\n";
#}#test

#check for exact kmer matches to all ref seqs
while(($kmerSeq,$kmerID) = each(%kmer)) {
#	while(($gb,$seq) = each(%ref)) {
	if (exists($kmer2{$kmerSeq})) { 
		#check for exact matches only
		#this can be modified later to allow Levenshtein distances (mismatches or gaps)
		#if necessary to account for sequence error or mutations
		$string = $kmer2{$kmerSeq};#$gb
		@string = split(',',$string);
		while($string[$i]) {
			$gb = $string[$i];
			if (exists($match{$gb})) {
				$string = $match{$gb};
				$newstring = $string.",".$kmerID;
				$match{$gb} = $newstring;
			}
			else {
				$match{$gb} = $kmerID;
			}
			$i++;
		}
		$i=0;
		$gb=();
		$string=();
		$newstring=();
#		print ".";
		$string=();
		@string=();
	}
#	print "\n";
}
print "Kmer matches found\n";

#for each reference sequence with more than 3 matches, count as present
while(($gb,$string) = each(%match)) {
	@kmerID = split(',',$string);
	$numKmer = scalar(@kmerID);
	if ($numKmer>=3) {#####edit this cutoff
		$present{$gb} = $ref{$gb};
		$count++;
	}
	@kmerID=();
	$numKmer=();
}
print "$count reference taxa were detected by at least 3 kmers (window=$windowSize, step=$stepSize)\n";

#print match table
open (OUT,">>","match.table") || die "Error cannot write to match.table: $!\n";
print OUT "GB of reference sequence\tkmerIDs that matched\n";

while(($gb,$string) = each(%match)){
	print OUT "$gb\t$string\n";
}
close OUT;

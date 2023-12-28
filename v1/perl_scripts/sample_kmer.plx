#!/usr/bin/perl
#March 27, 2012 by Terri Porter
#Script to grab kmers from cox gene sequences, need to input a window and step size from gb_seq.map (gb\tseq\n)
#usage perl sample_kmer.plx gb_seq.map

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb;
my $seq;
my $stepSize=50;#####edit this#####
my $windowSize=100;####edit this#####
my $freq;
my $length;
my $numKmer;
my $kmer;
my $j=0;

#declare array
my @map;
my @line;
my @seq;
my @kmer;
#my @error;
my @processed;

#declare hash
my %seq;
my %count;
my %processed;
my %error;

open (MAP,"<",$ARGV[0]) || die "Error cannot read gb_seq.map: $!\n";
@map = <MAP>;
close MAP;

#put gb and seq into seq hash, skip title line
$i=1;
while($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$gb = $line[0];
	$seq = $line[1];
	$seq{$gb} = $seq;
	$i++;
}
$i=0;

#go through seq hash, push all kmers into array
while(($gb,$seq) = each(%seq)) {
	if ($seq !~ /(B|[D-F]|[H-S]|[U-Z])/ig) {
#		$processed{$gb} = $seq;
		@seq = split("",$seq);
		$length = scalar(@seq);

		if ($length >= 500) { #####specify minimum reference sequence length here#####	
			$processed{$gb} = $seq;
			$numKmer = $length/$stepSize;
			$numKmer = int($numKmer);  #round down, only want full length kmers, no partials!
			while($i<$numKmer) {
				$kmer = substr($seq,$j,$windowSize);
				push(@kmer,$kmer);
				$i++;
				$j+=$stepSize;
			}
			$i=0;
			$j=0;
			@seq=();
			$length=();
			$numKmer=();
		}
		else {
			$error{$gb} = $seq;
		}
	}
	else {
		$error{$gb} = $seq;
	}
}

#use a hash to keep track of kmer frequency across all seqs
while ($kmer[$i]) {
	$kmer = $kmer[$i];
	if (exists($count{$kmer})) {
		$count{$kmer}++;
	}
	else {
		$count{$kmer}=1;
	}
	$i++;
}
$i=0;

#print file of kmer frequencies
open (OUT,">>","kmer.freq") || die "Error cannot write to kmer.freq: $!\n";

while (($kmer,$freq) = each(%count)) {
	print OUT "$kmer\t$freq\n";
}
close OUT;

#print all kmers to a fasta file for test assemblies
open (OUT2,">>","kmer.fasta") || die "Error cannot write to kmer.fasta: $!\n";

while($kmer[$i]) {
	$kmer = $kmer[$i];
	print OUT2 ">$i\n$kmer\n";
	$i++;
}
$i=0;
close OUT2;

#print list of gb seqs that contained ambiguities and were not processed
open (OUT3,">>","error.fasta") || die "Error cannot write to error.fasta: $!\n";

while (($gb,$seq) = each(%error)) {
	print OUT3 ">$gb\n$seq\n";
}
close OUT3;

#print list of gb seqs that were processed
open (OUT4,">>","processed.fasta") || die "Error cannot write to processed.fasta: $!\n";

while (($gb,$seq) = each(%processed)) {
	print OUT4 ">$gb\n$seq\n";
}
close OUT4;


#!/usr/bin/perl
#March 27, 2012 by Terri Porter
#April 30, 2012 edited to filter rbcL seqs, adding a taxonomy filter step to avoid sp.
#Script to grab kmers from cox gene sequences, need to input a window and step size from gb_seq.map (gb\tseq\n)
#NEW usage perl sample_kmer.plx gb_seq.map gb_org.map

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb;
my $org;
my $seq;
my $stepSize=5;#####edit this#####
my $windowSize=100;####edit this#####
my $lengthCutoff = 250;#####edit this#####
my $freq;
my $length;
my $numKmer;
my $kmer;
my $j=0;
my $scalar;

#declare array
my @map;
my @map2;
my @line;
my @seq;
my @kmer;
#my @error;
my @processed;
my @scalar;

#declare hash
my %seq;
my %org;
my %count;
my %processed;
my %error;

open (MAP,"<",$ARGV[0]) || die "Error cannot read gb_seq.map: $!\n";
@map = <MAP>;
close MAP;

#put gb and seq into seq hash
while($map[$i]) {
	$line = $map[$i];
	chomp $line;

	if ($line =~ /GB\tSequence/) { 
		$i++;
		next;
	}
	else {
		@line = split(/\t/,$line);
		$gb = $line[0];
		$seq = $line[1];
		$seq{$gb} = $seq;
	}
	
	$i++;
	$line=();
	@line=();
	$gb=();
	$seq=();
}
$i=0;

open (MAP2, "<", $ARGV[1]) || die "Error cannot read gb_org.map: $!\n";
@map2 = <MAP2>;
close MAP2;

#put gb and org into org hash, skip title line NEW
while ($map2[$i]) {
	$line = $map2[$i];
	chomp $line;

	if ($line =~ /GB\tOrganism/) {
		$i++;
		next;
	}
	else {
		@line = split(/\t/,$line);
		$gb = $line[0];
		$org = $line[1];
		$org{$gb} = $org;
	}

	$i++;
	$line=();
	@line=();
	$gb=();
	$org=();
}
$i=0;

#go through seq hash, push all kmers into array
while(($gb,$seq) = each(%seq)) {
	$org = $org{$gb};
	if ($org !~ /sp\./) {#avoid species without proper names, but keep hybrids
	
		if ($seq !~ /(B|[D-F]|[H-S]|[U-Z])/ig) {#avoid any ambiguities in sequence
#			$processed{$gb} = $seq;
			@seq = split("",$seq);
			$length = scalar(@seq);

			if ($length >= $lengthCutoff) { #####specify minimum reference sequence length here#####	
				$processed{$gb} = $seq;
				$numKmer = $length/$stepSize;
				$numKmer = int($numKmer);  #round down, only want full length kmers, no partials!
				while($i<$numKmer) {
					$kmer = substr($seq,$j,$windowSize);
					@scalar = split("",$kmer);
					$scalar = scalar(@scalar);
					if ($scalar == $windowSize) {#test to make sure kmer is full length
						push(@kmer,$kmer);
					}
					$i++;
					$j+=$stepSize;
					@scalar=();
					$scalar=();
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
open (OUT3,">>","bad_refs.fasta") || die "Error cannot write to bad_refs.fasta: $!\n";

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

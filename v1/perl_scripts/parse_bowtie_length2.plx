#!/usr/bin/perl
#Edited April 23, 2012 to require min. contiguous length instead of overlap
#Edited April 10, 2012 to require some kind of k-mer overlap in alignment for consensus
#March 30, 2012 by Terri Porter
#Script to parse bowtie -al output, create aligned fasta files without the reference, create consensus
#usage perl parse_bowtie.plx file.bowtie gb_seq.map

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign;

#declare var
my $i=0;
my $line;
my $gb;
my $start;
my $kmerSeq;
my $string;
my $newstring;
my $gbSeq;
my $length;
my $file;
my $startKmerSeq;
my $j=0;
my $kmerSeqLength;
my $partialLength;
my $gbLength;
my $lengthDifference;
my $str;
my $aln;
my $consensus;
my $numSeqs;
my $threshold;
my $lengthConsensus;
my $base;
my $count=0;
my $proportion;
my $k=0;
my $numParts;
my $part;
my $contiguousLength=50;#####edit cutoff here
my $minMapped = 1; #####edit minimum mapped here

#declare array
my @in;
my @line;
my @in2;
my @gbSeq;
my @string;
my @startKmerSeq;
my @kmerSeq;
my @output;
my @file;
my @consensus;
my @parts;
my @bases;

#declare hash
my %gb;#$gb=>$start|$kmerSeq csv
my %gbSeq;
my %gbLength;
my %part;

open (IN,"<",$ARGV[0]) || die "Error cannot read file.bowtie: $!\n";
@in = <IN>;
close IN;

#parse gb, start, and kmerString into hash for easier processing
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gb = $line[2];
	$start = $line[3];
	$kmerSeq = $line[4];
	
	if (exists($gb{$gb})) {
		$string = $gb{$gb};
		$newstring = $string.",".$start."|".$kmerSeq;
		$gb{$gb} = $newstring;
	}
	else {
		$gb{$gb} = $start."|".$kmerSeq;
	}

	$i++;
	@line=();
	$gb=();
	$start=();
	$kmerSeq=();
	$string=();
	$newstring=();
}
$i=0;

#parse gb_seq.map into hash, also make gb_length hash for easier searching
open (IN2,"<",$ARGV[1]) || die "Error cannot read gb_seq.map: $!\n";
@in2 = <IN2>;
close IN2;

$i=1;
while($in2[$i]) {#skip header line
	$line = $in2[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$gb = $line[0];
	$gbSeq = $line[1];
	$gbSeq{$gb} = $gbSeq;
	@gbSeq = split("",$gbSeq);
	$length = scalar(@gbSeq);
	$gbLength{$gb} = $length;

	$i++;
	$line=();
	@line=();
	$gb=();
	$gbSeq=();
	@gbSeq=();
	$length=();
}
$i=0;

#foreach gb, open a new aligned fasta file with aligned kmerStrings
while(($gb,$string) = each(%gb)) {
	$file = $gb.".fasta";
	
	open (OUT,">",$file) || die "Error cannot write to $file: $!\n";

	@string = split(',',$string);

	while($string[$i]) {
		$startKmerSeq = $string[$i];
		@startKmerSeq = split(/\|/,$startKmerSeq);
		$start = $startKmerSeq[0];
		$kmerSeq = $startKmerSeq[1];
		print OUT ">$i\n";

		while ($j<$start) {
			print OUT "-";
			$j++;
		}
		$j=0;

		print OUT "$kmerSeq";
		@kmerSeq = split("",$kmerSeq);
		$kmerSeqLength = scalar(@kmerSeq);
		$partialLength = $start+$kmerSeqLength;
		$gbLength = $gbLength{$gb};
		$lengthDifference = $gbLength-$partialLength;

		while ($lengthDifference>0) {
			print OUT "-";
			$lengthDifference--;
		}
		print OUT "\n";
		$i++;
		$startKmerSeq=();
		@startKmerSeq=();
		$start=();
		$kmerSeq=();
		@kmerSeq=();
		$kmerSeqLength=();
		$partialLength=();
		$gbLength=();
		$lengthDifference=();
	}
	$i=0;
	@string=();
	$file=();
	close OUT;
}

#for each gb.fasta file, create a consensus sequence and print to new file
@output = qx(ls | grep ".fasta");

open (OUT2,">>","consensus.seqs") || die "Error cannot write to consensus.seqs: $!\n";

while($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$gb = $file[0];

	#read in the alignment
	$str = Bio::AlignIO -> new('-file' => $file);
	$aln = $str->next_aln();#add it to hashes

	#grab conseus according to some threshold
	$numSeqs = $aln -> no_sequences;
	
	if ($numSeqs >= $minMapped) {#requires an alignment of more than one kmer!
		$consensus = $aln -> consensus_string(0);
		#instead of counting proportion of'?' bases, split on ? and count max. contiguous length instead
		@parts = split(/\?+/,$consensus);
		my $numParts = scalar(@parts);
		#print "Found $numParts parts in $gb\n";

		while ($parts[$j]) {
			my $part = $parts[$j];
			@bases = split("",$part);
			
			while($bases[$k]) {
				$base = $bases[$k];

				if ($base =~ /\w/) {
					$count++;
				}

				$k++;
				$base=();
			}
			$k=0;

			if ($count>0) {
				$part{$count} = $part;
			}

			$j++;
			$count=0;
			$part=();
			@bases=();
		}
		$j=0;
		@parts=();

#		while(my ($key,$value) = each (%part) ) {
#			print "$key => $value\n";
#		}
#		print "\n";

		for $count ( sort {$b<=>$a} keys %part) { #sort parts from largest to smallest
			if ($j == 0) {
				if ($count >= $contiguousLength) {
					$part = $part{$count};
					print OUT2 ">$gb\n$part\n";
				}
				$j++;
			}
		}
		$j=0;

	}

	$i++;
	$file=();
	@file=();
	$gb=();
	$str=();
	$aln=();
	$numSeqs=();
	$threshold=();
	$consensus=();
	%part=();
}
$i=0;
close OUT2;


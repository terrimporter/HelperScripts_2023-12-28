#!/usr/bin/perl
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
my $overlap=2;#####Edit this here, at least 2 kmers per column to be included in consensus
my $threshold;
my $lengthConsensus;
my $base;
my $count=0;
my $proportion;
my $ambiguityCutoff=50;#####Edit this here, allow up to 50 '?' columns

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

#declare hash
my %gb;#$gb=>$start|$kmerSeq csv
my %gbSeq;
my %gbLength;

open (IN,"<",$ARGV[0]) || die "Error cannot read file.bowtie: $!\n";
@in = <IN>;
close IN;

#parse gb, start, and kmerString into hash for easier processing
#while ($in[$i]) {
#	$line = $in[$i];
#	chomp $line;
#	@line = split(/\t/,$line);
#	$gb = $line[2];
#	$start = $line[3];
#	$kmerSeq = $line[4];
#	#print $kmerSeq."\n";#test works
#	if (exists($gb{$gb})) {
#		$string = $gb{$gb};
#		$newstring = $string.",".$start."|".$kmerSeq;
#		$gb{$gb} = $newstring;
#	}
#	else {
#		$gb{$gb} = $start."|".$kmerSeq;
#	}
#
	#print "@line\n";#test
	#print "gb:$gb\tstart:$start\tkmerSeq:$kmerSeq\n";
	#while(($gb,$string) = each(%gb)) {
	#	print "$gb\t=>\t$string\n";
	#}#test
#
#	$i++;
#	@line=();
#	$gb=();
#	$start=();
#	$kmerSeq=();
#	$string=();
#	$newstring=();
#}
#$i=0;

#parse gb_seq.map into hash, also make gb_length hash for easier searching
#open (IN2,"<",$ARGV[1]) || die "Error cannot read gb_seq.map: $!\n";
#@in2 = <IN2>;
#close IN2;

#$i=1;
#while($in2[$i]) {#skip header line
#	$line = $in2[$i];
#	chomp $line;

#	@line = split(/\t/,$line);
#	$gb = $line[0];
#	$gbSeq = $line[1];
#	print "gbSeq:$gbSeq\n";
#	$gbSeq{$gb} = $gbSeq;
#	@gbSeq = split("",$gbSeq);
#	$length = scalar(@gbSeq);
#	print "length:$length\n";
#	$gbLength{$gb} = $length;

#	$i++;
#	$line=();
#	@line=();
#	$gb=();
#	$gbSeq=();
#	@gbSeq=();
#	$length=();
#}
#$i=0;

#foreach gb, open a new aligned fasta file with aligned kmerStrings
#while(($gb,$string) = each(%gb)) {
#	$file = $gb.".fasta";
	
#	open (OUT,">",$file) || die "Error cannot write to $file: $!\n";

#	@string = split(',',$string);

#	while($string[$i]) {
#		$startKmerSeq = $string[$i];
		#print "starKmerSeq:$startKmerSeq\n";#test
#		@startKmerSeq = split(/\|/,$startKmerSeq);
#		$start = $startKmerSeq[0];
#		print "start:$start\n";
#		$kmerSeq = $startKmerSeq[1];
#		print "kmerSeq: $kmerSeq\n";
#		print OUT ">$i\n";

#		while ($j<$start) {
#			print OUT "-";
#			$j++;
#		}
#		$j=0;

#		print OUT "$kmerSeq";
#		@kmerSeq = split("",$kmerSeq);
#		$kmerSeqLength = scalar(@kmerSeq);
#		$partialLength = $start+$kmerSeqLength;
#		$gbLength = $gbLength{$gb};
#		$lengthDifference = $gbLength-$partialLength;

#		while ($lengthDifference>0) {
#			print OUT "-";
#			$lengthDifference--;
#		}
#		print OUT "\n";
#		$i++;
#		$startKmerSeq=();
#		@startKmerSeq=();
#		$start=();
#		$kmerSeq=();
#		@kmerSeq=();
#		$kmerSeqLength=();
#		$partialLength=();
#		$gbLength=();
#		$lengthDifference=();
#	}
#	$i=0;
#	@string=();
#	$file=();
#	close OUT;
#}

#for each gb.fasta file, create a consensus sequence and print to new file

@output = qx(ls | grep ".fasta");

open OUT2,">>consensus_min2_max50.seqs.quick" || die "Error cannot write to consensus.seqs: $!\n";

while($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$gb = $file[0];
#	print "$file\n";

	#read in the alignment
	$str = Bio::AlignIO -> new('-file' => $file);
	$aln = $str->next_aln();#add it to hashes

	#grab conseus according to some threshold
	$numSeqs = $aln -> no_sequences;
	#print "numSeqs: $numSeqs\n";#test
	#$threshold = ($overlap/$numSeqs)*100;
	
	if ($numSeqs >= $overlap) {#requires an alignment of more than one kmer!
		$threshold = (($overlap/$numSeqs)*100);#do not use int
		#default = char needs to be present at least 0% to appear, ambiguity codes used, lowercase for DNA using AlignIO, then lowercase if gaps in alignment using SimpleAlign
	#	print "threshold: $threshold\n";#test
		$consensus = $aln -> consensus_string($threshold);
	#	print "consensus:$consensus\n";#test
		@consensus = split("",$consensus);
		$lengthConsensus = scalar(@consensus);
		while ($consensus[$j]) {
			$base = $consensus[$j];
			
			if ($base eq "?") {
				$count++;
			}
			$j++;
		}
		$j=0;
		$proportion = ($count/$lengthConsensus)*100;#do not use int
		#print "proportionGap:$proportion\n";#test

		$count=0;#RESET!!!!!!!
	#	print "numGap:$count\n";
	#	print "length:$lengthConsensus\n";
	#	print "proportionGap:$proportion\n";#test
		if ($proportion <= $ambiguityCutoff) { #only print to file if less than 50% ambituous bases
			print OUT2 ">$gb\n$consensus\n";
	#		print OUT2 "test\n";
			print ".";
			$proportion=();
		}
			
		#default = char needs to be present at least 0% to appear, ambiguity codes used, lowercase for DNA using AlignIO, then lowercase if gaps in alignment using SimpleAlign
	}
	print "\n";
	$i++;
	$file=();
	@file=();
	$gb=();
	$str=();
	$aln=();
	$numSeqs=();
	$threshold=();
	$consensus=();
	@consensus=();
	$lengthConsensus=();
	$base=();
}
$i=0;
close OUT2;

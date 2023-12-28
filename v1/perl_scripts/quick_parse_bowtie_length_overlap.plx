#!/usr/bin/perl
#Edited May 7, 2012 to parse fasta files already available
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
my $file;
my $gb;
my $j=0;
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
my $overlap=2;#####edit overlap here

#declare array
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

#for each gb.fasta file, create a consensus sequence and print to new file
@output = qx(ls | grep ".fasta");

open (OUT2,">>","consensus.50contig.min2ol.seqs.quick") || die "Error cannot write to outfile: $!\n";

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
	
	if ($numSeqs >= 1) { 
		#requires an alignment of more than one kmer!
		$threshold = ($overlap/$numSeqs) * 100; 
		#no int
		$consensus = $aln -> consensus_string($threshold); 
		#still require certain amount of overlap!
		#instead of counting proportion of'?' bases, split on ? and count max. contiguous length instead
		@parts = split(/\?+/,$consensus);
		#print "arrayParts:@parts\n";
		$numParts = scalar(@parts);
		#print "Found $numParts parts in $gb\n";
	
		if ($numParts == 0 ) {
			#print "no parts to parse\n";
		}
		else {
			print "try to parse parts\n";
			
			foreach $part (@parts) {
				print "part:$part\n";
				@bases = split("",$part);
				
				if (@bases) {
					while($bases[$k]) {
						$base = $bases[$k];
	
						if ($base =~ /\w/) {
							$count++;
							#print "$base";
						}
	
						$k++;
						$base=();
					}
				}
				@bases=();
				#print "\n";
				$k=0;

				if ($count>0) {
					$part{$count} = $part;
				}

				$j++;
				$count=0;
				$part=();
			}
			$j=0;
			@parts=();
		}
		$numParts=();
		
		if (%part) {
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
		%part=();
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
}
$i=0;
close OUT2;


#!/usr/bin/perl

#Edited Aug. 2, 2012 for Bendan, using aligned muscle fasta files
#June 6, 2012 by Terri Porter
#Script to create a consensus sequence from an aligned fasta file
#usage perl create_consensus.plx

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign;

#declare var
my $i=0;
my $file;
my $basename;
my $str;
my $aln;
my $numSeqs;
my $minSeqs = 2;
my $consensus;
my $threshold = 0;### BRENDAN EDIT THIS IF YOU WANT ### 

#consensus threshold can range from 0 to 100 (%). BIO::SimpleAlign default is 0.  
#consensus residue has to appear in at least x percent of the aligned sequences otherwise '?' will be put in consensus.

#declare array
my @output;
my @file;

@output = qx(ls | grep 'out\$'); #only work with filenames that end with .out

open (OUT, ">>", "consensus.seqs") || die "Error cannot open outfile: $!\n";

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	@file = split(/\./, $file);
	$basename = $file[0];

	#read in the alignment
	$str = Bio::AlignIO -> new ('-file' => $file);
	$aln = $str -> next_aln();

	#grab consensus according to some threshold
	$numSeqs = $aln -> no_sequences;
	if ($numSeqs >= $minSeqs) {

		$consensus = $aln -> consensus_string($threshold);
	
		if ($consensus =~ /^[?]+/) { #remove starting ?'s
			$consensus =~ s/^[?]+//;
		}
		
		if ($consensus =~ /[?]+$/) { #remove trailing ?'s
			$consensus =~ s/[?]+$//;
		}


		if ($consensus =~ /\?/) { #convert internal ?'s to N's
			$consensus =~ s/\?/N/g;
		}
		print OUT ">$basename\n$consensus\n";
	}
	$i++;
	$file=();
	@file=();
	$basename=();
	$str=();
	$aln=();
	$numSeqs=();
	$consensus=();
}
$i=0;
close OUT;

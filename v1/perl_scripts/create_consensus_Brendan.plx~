#!/usr/bin/perl
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
my $gb;
my $str;
my $aln;
my $numSeqs;
my $minMapped = 1; ### Edit number of mapped baits here
my $consensus;

#declare array
my @output;
my @file;

@output = qx(ls | grep 'sim\$');

open (OUT, ">>", "consensus.seqs") || die "Error cannot open outfile: $!\n";

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	@file = split(/\./, $file);
	$gb = $file[0];
	#print "$gb\n";
	#read in the alignment
	$str = Bio::AlignIO -> new ('-file' => $file);
	$aln = $str -> next_aln(); #add it to hashes

	#grab consensus according to some threshold
	$numSeqs = $aln -> no_sequences; #get number of sequences in alignment
	if ($numSeqs >= $minMapped) {

		$consensus = $aln -> consensus_string();
		#print "$consensus\n";
		if ($consensus =~ /^[?]+/) { #remove starting ?'s
			$consensus =~ s/^[?]+//;
			#print "start $consensus\n";
		}
		
		if ($consensus =~ /[?]+$/) { #remove trailing ?'s
			$consensus =~ s/[?]+$//;
			#print "end $consensus\n";
		}


		if ($consensus =~ /\?/) { #convert internal ?'s to N's
			$consensus =~ s/\?/N/g;
			#print "internal $consensus\n";
		}
		print OUT ">$gb\n$consensus\n";
		#print "\n";
	}
	$i++;
	$file=();
	@file=();
	$gb=();
	$str=();
	$aln=();
	$numSeqs=();
	$consensus=();
}
$i=0;
close OUT;

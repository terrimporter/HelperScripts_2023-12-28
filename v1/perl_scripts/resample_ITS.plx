#!/usr/bin/perl
#Feb.7,2011 by Terri Porter
#Script to create new fragmented ITS datasets by grabbing short fragments of the 5' and 3' ends in increasing lengths starting at 50bp going up to 600bp if possible
#Use these new datasets for blast and MEGAN, and SAP analysis then compare results with results from the complete ITS dataset
# usage $perl $resample_ITS.plx ITS.fasta

use strict;
use warnings;

#var
my $i=0;
my $fragment_size;
my $outfile_name;
my $j=0;
my $line;
my $header;
my $sequence;
my $k=0;
my $base;
my $cutoff;
my $sequence_fragment;

#array
my @fasta;
my @fragment_sizes;
my @sequence;
my @sequence_fragment;
my @sequence_reverse;
my @sequence_fragment_reverse;

open (IN,"<",$ARGV[0]) || die ("Error cannot read from ITS.fasta: $!\n");
@fasta = <IN>;
close IN;
foreach $line (@fasta) {
	chomp $line;
}
#print "fasta array:\n@fasta\n";#test

@fragment_sizes = (50,100,150,200,250,300,350,400,450,500,550,600);
#print "fragment sizes array:\n@fragment_sizes\n";#test

#create 5' fragments

while ($fragment_sizes[$i]) {
	$fragment_size = $fragment_sizes[$i];

	$outfile_name = "5_prime_ITS_".$fragment_size.".fasta";
	open (OUT,">>",$outfile_name) || die ("Error cannot write to outfile: $!\n");
	
	while ($fasta[$j]) {
		$line = $fasta[$j];
		
		if ($line =~ /^>/) {
			#print "found header\n";#test
			$header = $line;
			print OUT "$header\n";
		}
		else {
			#print "found sequence\n";#test
			$sequence = $line;
			@sequence = split(//,$sequence);
			#print "sequence array:\n@sequence\n";#test
			
			while ($sequence[$k]) {
				$base = $sequence[$k];
				$cutoff = $fragment_size-1;				
				if ($k<=$cutoff) {
					push(@sequence_fragment, $base);
				}
				$k++;
			}
			$k=0;
			$sequence_fragment = join("",@sequence_fragment);
			print OUT "$sequence_fragment\n";			
			@sequence_fragment=();
		}
		$j++;
	}
	close OUT;
	$j=0;
	$i++;
}

$i=0;

while ($fragment_sizes[$i]) {
	$fragment_size = $fragment_sizes[$i];

	$outfile_name = "3_prime_ITS_".$fragment_size.".fasta";
	open (OUT,">>",$outfile_name) || die ("Error cannot write to outfile2: $!\n");

	while ($fasta[$j]) {
		$line = $fasta[$j];

		if ($line =~ /^>/) {
			$header = $line;
			print OUT "$header\n";
		}
		else {
			$sequence = $line;
			@sequence = split(//,$sequence);
			@sequence_reverse = reverse(@sequence);

			while ($sequence_reverse[$k]) {
				$base = $sequence_reverse[$k];
				$cutoff = $fragment_size-1;
				if ($k<=$cutoff) {
					push(@sequence_fragment_reverse, $base);
				}
				$k++;
			}
			$k=0;
			@sequence_fragment = reverse(@sequence_fragment_reverse);
			$sequence_fragment = join("",@sequence_fragment);
			print OUT "$sequence_fragment\n";
			@sequence_fragment_reverse=();
			@sequence_fragment=();
		}
		$j++;
	}
	close OUT;
	$j=0;
	$i++;
}

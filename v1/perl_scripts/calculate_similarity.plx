#!/usr/bin/perl
#Nov.19,2010 by Terri Porter
#Script to calculate average % similarity for an aligned region
#Requires aligned primers and aligned seqs of the same length
#usage $perl calculate_similarity.plx name.primers2 name.seqs

use strict;
use warnings;

#declare var
my $i=0;
my $current_seq;
my $length;
my $j=0;
my $current_primer;
my $k=0;
my $compare_seq;
my $compare_primer;
my $l=0;
my $similarity;
my $m=0;
my $line;
my $num_elements;
my $sum_similarities;
my $ave_similarity;

#declare array
my @primers;
my @seqs;
my @current_sequence;
my @current_primers;
my @temp;
my @line;

open (PRIMER,"<",$ARGV[0]) || die ("Error: $!\n");
@primers = <PRIMER>;
close PRIMER;

open (SEQ,"<",$ARGV[1]) || die ("Error: $!\n");
@seqs = <SEQ>;
close SEQ;

open (TEMP,">>","similarity.temp") || die ("Error: $!\n");

while ($seqs[$i]) {
	$current_seq = $seqs[$i];
	chomp $current_seq;
	@current_sequence = split(//,$current_seq);
	$length = scalar(@current_sequence);

	while ($primers[$j]) {
		$current_primer = $primers[$j];
		chomp $current_primer;
		@current_primers = split(//,$current_primer);
		
		while ($current_sequence[$k]) {
			$compare_seq = $current_sequence[$k];
			$compare_primer = $current_primers[$k];
			if ($compare_seq eq $compare_primer) {
				$l++;
			}
			$k++;
		}
		$similarity = $l/$length*100;
		print TEMP "$similarity\t";
		$j++;
		$l=0;
		$k=0;
	}
	print TEMP "\n";
	$j=0;
	$i++;
}
close TEMP;

open (TEMP,"<","similarity.temp") || die ("Error: $!\n");
@temp = <TEMP>;
close TEMP;

open (OUT,">>","ave_similarity.txt") || die ("Error: $!\n");

while ($temp[$m]) {
	$line = $temp[$m];
	chomp $line;
	@line = split(/\t/,$line);
	$num_elements = scalar(@line);
	$sum_similarities =0;
	($sum_similarities+=$_) for @line;
	$ave_similarity = $sum_similarities/$num_elements;
	print OUT "$ave_similarity\n";
	$m++;
}

unlink("similarity.temp");

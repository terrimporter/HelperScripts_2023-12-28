#!/usr/bin/perl
#January 4, 2013 edit to skip blank lines in some of Joel's fasta files
#October 10, 2009
#Script to parse a fasta file to obtain sequence lengths for each sequence.  This can be opened in excel to calculate mean, mode, range, etc.
#USAGE $perl fasta_stats.plx < infile.txt >outfile.txt
#modified Nov.17, 2010 to count number of records in file

use strict;
use warnings;

use lib '/home/tp45/modules';
use Statistics::Lite qw(:all);

#declare variables
my $line;
my $flag=0;
my $seq1;
my $seq2;
my $seq;
my $i=0;
my $length;
my $min;
my $max;
my $mean;
my $mode;

#declare array
my @seq;
my @split;
my @length;

while(<>){
	$line = $_;
	chomp($line);

	if ($flag==0){

		if (/>/){
			next;
		}
		elsif (/^s+/) { ### skip blank lines ###
			next;
		}
		else {
			$seq1 = $line;
			$flag=1;
		}
		
	}	
	elsif ($flag==1){
		
		if (/>/){
			$flag=0;
			push (@seq, $seq1);
		}
		elsif (/^s+/) { ### skip blank lines ###
			next;
		}
		else {
			$seq2 = $line;
			$seq1 = $seq1.$seq2;
		}
	}
}
push (@seq, $seq1);#don't forget to add last seq in file!
	
#test-passed
#foreach (@seq){
#	my $element = $_;
#	print $element,"\n";
#}

while ($seq[$i]){
	@split = split(//, $seq[$i]);
	$length = scalar(@split);
	push (@length, $length);
	@split =();#empty array
	$i++;
}

#test-passed
foreach (@length){
	my $element = $_;
	print $element, "\n";
}

$min = min (@length);
$max = max (@length);
$mean = mean (@length);
$mode = mode (@length);
my $num = scalar(@seq);

print "NumSeqs\tMin\tMax\tMean\tMode\n";
print $num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";

#!/usr/bin/perl
#Terri Porter, October 6, 2017
#Script to grab seq stats from fastq files
#USAGE perl fastq_stats.plx infile.fastq.gz

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $num;
my $i=0;
my $line;
my $length;
my $num2;
my $count;
my $min;
my $max;
my $mean;
my $median;
my $mode;

#declare array
my @allseqs;
my @line;
my @lengths;

@allseqs = `zcat $ARGV[0] | awk '(NR%4==2)'`;
$num = scalar (@allseqs);

while ($allseqs[$i]) {
	$line = $allseqs[$i];
	chomp $line;

	@line = split(//,$line);
	$length = scalar @line;
	push @lengths, $length;
	$i++;
}
$i=0;

$num2 = scalar (@lengths);

if ($num != $num2) {
	print "Possible error\n";
}

$count = count (@lengths);
$min = min (@lengths);
$max = max (@lengths);
$mean = mean (@lengths);
$median = median (@lengths);
$mode = mode (@lengths);

print STDOUT $ARGV[0]."\t".$count."\t".$min."\t".$max."\t".$mean."\t".$median."\t".$mode."\n";

$i=0;

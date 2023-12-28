#!/usr/bin/perl
#Terri Porter, Aug. 20, 2010
#Script to count average number of characters in temporary phylip file
#Usage $perl count_phylip_char.plx phylip_temp.txt

use strict;
use warnings;

#declare variables
my $line;
my $char;
my $j=0; #counter
my $k=0; #counter
my $sum;
my $average;

#declare array
my @chars;
my @num_char;

open (IN, '<', 'phylip_temp.txt') || die ("Error: Cannot open phylip_temp file $!");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /\w+/) {
		$line =~ /\d+\s+(\w+)/;
		$char = $1;
		@chars = split (//, $char);
		foreach (@chars) {
			$k++;
		}
		push (@num_char, $k);
		$k=0;
		$j++;
	}
}

foreach (@num_char) {
	$sum += $_;
}

$average = $sum/$j;
print $average;

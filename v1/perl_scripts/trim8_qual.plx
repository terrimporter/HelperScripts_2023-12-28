#!/usr/bin/perl
#Written by Terri Porter, latest update April 2, 2013, for Porter et al., 2013 MER
#Script to remove phred scores associated with MID (10 bp) in file.qual
#USAGE $perl trim10_qual.plx < file.qual > file.qual.MIDremoved

use strict;
use warnings;

#declare var
my $line;
my $header;
my $seq;
my $i=1;
my $trimmed_seq;

#declare array
my @line;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
	}
	else {
		$seq = $line;
		@line = split (/ /,$line);##### split by space instead #####
		while ($i <= 9) { ### Nagissa / MR DNA used 8 bp barcode tags ###
			shift(@line);
			$i++;
		}
		$trimmed_seq = join(' ',@line); ##### join by space instead #####
		print "$header\n$trimmed_seq\n";	
		$i=1;
	}
}	

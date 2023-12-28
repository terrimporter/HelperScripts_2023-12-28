#!/usr/bin/perl
#Jan.3, 2012 edited to remvoe first 10 phred scores from *.qual instead
#Nov.16, 2010 by Terri Porter
#Script to remove first 10 base pairs from fasta file sorted by MID (10bp)
#usage $perl trim10.plx < fasta.file > trimmed_fasta.file

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
		@line = split (/ /,$line);##### split by pace instead #####
		while ($i <= 10) {
			shift(@line);
			$i++;
		}
		$trimmed_seq = join(' ',@line); ##### join by space instead #####
		print "$header\n$trimmed_seq\n";	
		$i=1;
	}
}	

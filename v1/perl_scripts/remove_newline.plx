#!/usr/bin/perl
#Terri Porter, edited August 3, 2016 to work right
#Terri Porter, August 19, 2010
#Script to remove newline character at the end of a sequence in a fasta formatted file
#Usage $perl remove_newline.plx < infile > outfile

use strict;
use warnings;

#declare variables
my $line;
my $i=0; #flag first line

while (<>) {
	$line = $_;
	chomp $line;

	if ($i==0) { #first line
		if ($line =~ /^>/) {
			print $line,"\n";
			$i++;
		}
	}
	elsif ($i > 0) {

		if ($line =~ /^>/) {
			print "\n",$line,"\n";
		}
		elsif (eof()) { #check for last line
			print $line,"\n";
		}
		else {
			print $line;
		}
	}
}

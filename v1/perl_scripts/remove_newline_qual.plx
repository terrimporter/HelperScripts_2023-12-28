#!/usr/bin/perl
#Terri Porter, August 19, 2010
#Script to remove newline character at the end of a sequence in a fasta formatted file
#Usage $perl remove_newline.plx < infile > outfile

use strict;
use warnings;

#declare variables
my $line;
my $flag=0;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		print "\n",$line,"\n";
		$flag=0;
	}
	elsif ($flag==0) {
		print $line;
		$flag=1;
	}
	elsif ($flag>0) {
		print " $line";
	}
}

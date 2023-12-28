#!/usr/bin/perl

#Terri Porter, August 19, 2010

#Script to grab just COI sequences from the fasta formatted Dasmahapatra dataset from gengank
#Usage $perl grab_COI_only.plx < infile > outfile

use strict;
use warnings;

#declare some variables
my $line;
my $fastaheader=0; #flag
my $COIsequence=0; #flag

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$fastaheader = 1;
		if ($line =~ /(COI)/) {
			$COIsequence = 1;
			if ($fastaheader eq 1 && $COIsequence eq 1) {
				print "\n",$line,"\n";
				next;
			}
		}
		if ($line !~ /(COI)/) {
			$COIsequence = 0;
			next;
		}		
	}
	elsif ($line !~ /^>/) {
		$fastaheader =0;
		if ($COIsequence eq 1) {
			print $line;
			next;
		}
		elsif ($COIsequence eq 0) {
			next;
		}
	}
}


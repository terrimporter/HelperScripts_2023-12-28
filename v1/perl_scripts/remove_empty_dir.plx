#!/usr/bin/perl
#Terri Porter, Aug. 22, 2010
#Script to remove empty sequences from the seqtrim outfile, then run fasta_stats to get stats after trimming
#Usage $perl remove_empty.plx < in > out

use strict;
use warnings;

#declare variables
my $line;

#declare array
my @line;

my $dir = "/home/terri/Mehrdad_soil_grant/second_batch/subsamples/from_seqtrim/";
opendir DH, $dir;
my @files = readdir (DH);

foreach my $files (@files) {
	if ($files eq "." ) {
		next;
	}
	elsif ($files eq "..") {
		next;
	}
	else {
		my $path_to_infile = $dir.$files;
		open (IN, "<", $path_to_infile) || die ("Error:$!\n");

		while (<IN>) {
			$line = $_;
			chomp $line;
			if ($line =~ /^>/) {
				push (@line, $line); #add to end of array
			}
			else {
				if ($line =~ /\w+/) {
					push (@line, $line);
				}
				else {
					pop (@line);
				}
			}
		}
		my $outfile=$files.".noempty";
		open (OUT,">>",$outfile) || die ("Error:$!\n");
		foreach (@line) {
			$line = $_;
			print OUT $line."\n";
		}
		close IN;
		close OUT;
		@line=();
	}
}

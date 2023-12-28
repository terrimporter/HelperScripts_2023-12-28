#!/usr/bin/perl
#Oct.13,2010 by Terri Porter
#Script to reformat fasta headers after sorting by abundance
#usage $perl reformat_cluster_header.plx < in > out
#Dec.6,2010 edited to work properly.
#March 3,2011 edited to print only readID in header
#
use strict;
use warnings;

#declare var
my $line;
my $readID;
my $symbol;
my $rest;
my $readIDpart;

#declare array
my @line;
my @readIDpart;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		if ($line =~ />/) {
			$line =~ s/>//;
			@line = split(/ /,$line);
			$readIDpart = $line[0];
			@readIDpart = split(/\|/,$readIDpart);
			$readID = $readIDpart[2];
			print ">$readID\n";
		}
	}
	else {
		print "$line\n";
	}
}

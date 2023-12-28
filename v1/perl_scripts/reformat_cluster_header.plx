#!/usr/bin/perl
#Oct.13,2010 by Terri Porter
#Script to reformat fasta headers after sorting by abundance
#usage $perl reformat_cluster_header.plx < in > out
#Dec.6,2010 edited to work properly.
#
use strict;
use warnings;

#declare var
my $line;
my $clusterID;
my $symbol;
my $rest;

#declare array
my @line;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		if ($line =~ />\d+\|\*\|\w{14}/) {
			@line = split(/\|/,$line);
			$clusterID = $line[2];
			$symbol = $line[3];
			$rest = $line[4];
			print ">$clusterID|$symbol|$rest\n";
		}
	}
	else {
		print "$line\n";
	}
}

#!/usr/bin/perl
#Terri Porter, Aug.23, 2010
#Script to add qutoes around taxon names for bpp.phy files so they open properly in Mesquite
#Usage $perl add_quotes.plx < in > out

use strict;
use warnings;

#declare variables
my $line;
my $label;
my $seq;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /\d+\t\d+/) {
		$line =~ /(\d+\t\d+)/;
		$line = $1;
		print "\t$line\n";
		next;
	}
	elsif ($line =~ /\d{9}.{2,3}\t\D+/) {
		$line =~ /(\d{9}.{2,3})\t(\D+)/;
		$label = $1;
		$seq = $2;
		print "'$label'\t$seq\n";
		next;
	}

}

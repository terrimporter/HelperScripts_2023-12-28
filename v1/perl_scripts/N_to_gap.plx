#!/usr/bin/perl
#Terri Porter, Aug.20, 2010
#Script to change N's to gaps
#Usage: perl N_to_gap.plx <in > out

use strict;
use warnings;

#declare var
my $line;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /N/) {
		$line =~ s/N/-/g;
		print "$line\n";
	}
	else {
		print "$line\n";
	}
}

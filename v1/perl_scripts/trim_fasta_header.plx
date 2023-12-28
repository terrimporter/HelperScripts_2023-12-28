#!/usr/bin/perl
#Sept. 20, 2010 by Terri Porter
#Script to make fasta headers more concise for SAP processing
#usage $perl trim_fasta_header.plx < in > out

use strict;
use warnings;

#declare variables
my $line;
my $id;

#declare arrays
my @line;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ />/) {
		@line = split (/\|/, $line);
		$id = $line[0];
		print "$id\n";
	}
	else {
		print "$line\n";
	}
}

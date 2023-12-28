#!/usr/bin/perl
#Written by Terri Porter, lastest updat April 2, 2013, for Porter et al., 2013 MER
#Script to reformat fasta header to keep readid only
#USAGE $perl reformat_fsata_id.plx < infile > outfile

use strict;
use warnings;

#var
my $i=0;
my $line;
my $header;

#arrays
my @in;

@in = <STDIN>;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line  =~ /^(>\S+)\s+/;
		$header = $1;
		print STDOUT $header."\n";
	}
	else {
		print STDOUT $line."\n";
	}
	$i++;
	$line=();
	$header=();
}
$i=0;

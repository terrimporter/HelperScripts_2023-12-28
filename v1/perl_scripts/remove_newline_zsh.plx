#!/usr/bin/perl
#Sept. 2, 2016, edit to work with run_removenewline.sh on a directory of fasta files
#Terri Porter, edited August 3, 2016 to work right
#Terri Porter, August 19, 2010
#Script to remove newline character at the end of a sequence in a fasta formatted file
#Usage $perl remove_newline.plx infile outfile

use strict;
use warnings;

#declare variables
my $line;
my $i=0; #flag first line
my $j=0;

#declare array
my @fasta;

open (FASTA, "<", $ARGV[0]) || die "Error cannot open file: $!\n";
@fasta = <FASTA>;
close FASTA;

open (OUT, ">>", $ARGV[1]) || die "Error cannot open outfile: $!\n";

while ($fasta[$j]) {
	$line = $fasta[$j];
	chomp $line;

	if ($i==0) { #first line
		if ($line =~ /^>/) {
			print OUT $line,"\n";
			$i++;
			$j++;
		}
		else {
			$j++;
		}
	}
	elsif ($i > 0) {

		if ($line =~ /^>/) {
			print OUT "\n",$line,"\n";
			$j++;
		}
		elsif (eof()) { #check for last line
			print OUT $line,"\n";
			$j++;
		}
		else {
			print OUT $line;
			$j++;
		}
	}
}
$j=0;
close OUT;

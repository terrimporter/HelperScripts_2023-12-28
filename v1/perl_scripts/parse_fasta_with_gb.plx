#!/usr/bin/perl
#April 12, 2011 by Terri Porter
#Script to parse ITS.fasta using a list of gb accessions (no version)
#usage $perl parse_fasta_with_gb.plx gb.list ITS.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $j=0;
my $gb;
my $next_line;

#declare array
my @gb;
my @fasta;

open (LIST,"<",$ARGV[0]) || die ("Error cannot read list: $!\n");
@gb = <LIST>;
close LIST;

open (FASTA,"<",$ARGV[1]) || die ("Error cannot read fasta: $!\n");
@fasta = <FASTA>;
close FASTA;

open (OUT,">>","fasta") || die ("Error cannot write to fasta: $!\n");

while ($gb[$j]) {
	$gb = $gb[$j];
	chomp $gb;

	while ($fasta[$i]) {
		$line = $fasta[$i];
		chomp $line;
	
		if ($line =~ /$gb/) {
			print OUT "$line\n";
			$i++;
			$next_line = $fasta[$i];
			chomp $next_line;
			print OUT "$next_line\n";
			last;
		}
		$i++;
	}
	$i=0;
	$j++;
}
close OUT;

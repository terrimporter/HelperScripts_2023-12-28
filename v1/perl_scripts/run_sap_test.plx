#!/usr/bin/perl
#March 15,2011 by Terri Porter
#Script to use a mapped fasta file (with gi numbers in the header) to conduct sap searches using the forceexcludegilist option.
#usage $perl run_sap_test.plx ITS.fasta.mapped

#declare var
my $line;
my $i=0;
my $gi;
my $filename;
my $j=0;

#declare array
my @fasta;
my @line;
my @gi;
my @output;

use strict;
use warnings;

open (FASTA,"<",$ARGV[0]) || die ("Error cannot read from mapped fasta file: $!\n");
@fasta = <FASTA>;
close FASTA;

print "\nCreating individual gi.fasta files\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ s/^>//;
		@line = split(/\|/,$line);
		$gi = $line[0];
		push(@gi,$gi);
		$filename = $gi.".fasta";
	
		open (OUT,">>",$filename) || die ("Error cannot write to fasta outfile: $!\n");
		print OUT ">$line\n";
	}
	else {
		print OUT "$line\n";
		close OUT;
	}
	$i++;
}

print "\nStarting sap runs\n";

while ($gi[$j]) {
	$gi = $gi[$j];
	$filename = $gi.".fasta";####edit settings below###
	system("sap $filename --project NJrun --assignment Barcoder --minidentity 0.70 --phyla 1 --classes 1 --orders 1 --families 1 --genera 3 --individuals 5 --subspecieslevel --forceexcludegilist '$gi' ");
	$j++;
}

#!/usr/bin/perl
#Oct. 31, 2016 by Terri Porter
#Script to remove bacterial sequences from testNBC.fasta before screening for bacterial contaminants using BLAST
#USAGE perl filter_Bact_from_fasta.plx testNBC.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;

#declare array
my @in;

open (IN, "<", $ARGV[0]) || die "Cannot open infile:$!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "testNBC_noBact.fasta") || die "Cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		if ($line =~ /Bacteria/) {
			$i = $i+2;
			next;
		}
		else {
			print OUT $line."\n";
		}

	}
	else {
		print OUT $line."\n";
	}
	$i++;
}
$i=0;

close OUT;

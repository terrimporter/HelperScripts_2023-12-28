#!/usr/bin/perl
#Nov.12,2010 by Terri Porter
#Scrip to turn agrep output into fasta formatted file
#usage $perl agrep_to_fasta.plx agrep.out

use strict;
use warnings;

#declare var
my $line;
my $seq;
my $id;

#declare array
my @line;

open (IN,"<",$ARGV[0]) || die "Error: $!\n";

open (OUT,">","agrep.fasta") || die "Error: $!\n";

while (<IN>) {
	$line = $_;
	chomp $line;
	@line = split(/\|/,$line);
	$seq = $line[0];
	$id = $line[1];
	print OUT ">$id\n$seq\n";
}

close IN;
close OUT;

#!/usr/bin/perl
#Oct.6,2010 by Terri Porter
#Script to filter clone_summary.txt to keep only the best taxon assignment with highest posterior probability
#usage $perl filter_clone_summary.txt infile

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $accession;
my $genus;
my $species;
my $pp;
my $j;
my $best_pp;
my $x;
my $index;
my $newline;

#declare array
my @line;
my @file;
my @accession;
my @genus;
my @species;
my @pp;
my @pp2;
my @newline;

open (IN, "<",$ARGV[0]) || die ("Error:$!\n");
@file = <IN>;
foreach $x (@file) {
	chomp $x;
}
close IN;

while ($file[$i]) {
	$line = $file[$i];
	@line = split(/\t/,$line);
	$accession = $line[0];
	push (@accession,$accession);
	$genus = $line[1];
	$species = $line[2];
	$pp = $line[3];
	if ($pp =~ /~/) {
		$pp =~ s/~//;
	}
	$newline = "$genus\t$species\t$pp";
	push(@newline,$newline);
	push (@pp,$pp);
	$j=$i+1;
	if ($file[$j] =~ /$accession/) {
		$i++;
		next;
	}
	else {
		@pp2 = sort {$a <=> $b} @pp;
		$best_pp = pop(@pp2);
		foreach $newline (@newline) {
			if ($newline =~ /$best_pp/) {
				print "$accession\t$newline\n";
			}
		}
		@newline=();
		@accession=();
		@genus =();
		@species=();
		@pp=();
		@pp2=();
		$i++;
	}
}	

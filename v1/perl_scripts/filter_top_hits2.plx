#!/usr/bin/perl
#Script edited so it grabs the last good hit and prints it!!! use this one instead of filter_top_hits.plx
#Modified Oct.15,2010 to sort by QueryID, sort by bitscore, and keep only top hit if bitscore > 100 (ca. 80% of smallest target region)
#Oct.14,2010 by Terri Porter
#Script to parse top blast hits and keep entries with 99-100% sequence identity, corresponding to up to 1 mismatch or 1 gap opening allowed
#usage $perl filter_top_hits.plx blastall.out > filtered_blastall.out

use warnings;
use strict;

#declare variables
my $i=0;
my $line;
my $identity;
my $queryID_current;
my $queryID_previous = "nil";
my $flag=0;
my $value;
my $key;
my $j=0;
my $bitscore_entry;
my $values_index;
my $line_index;
my $best_line;
my $best_bitscore;
my $x;
my $size;

#declare array
my @file;
my @line;
my @temp;
my @values;
my @keys;
my @entry;

#declare hashes
my %bitscore;
my %line;

open (IN,"<",$ARGV[0]) || die ("Error: $!\n");
@file = <IN>;
close IN;

#filter for %identities >= 99%
while ($file[$i]) {
	$line = $file[$i];
	chomp $line;
	@line =split(/\t/,$line);
	$identity = $line[2];
	if ($identity >= 99) {######modify percent identity cutoff here########
		push(@temp,$line);
	}
	$i++;
}

#for each query, filter by bitscore
$i=0;
while ($temp[$i]) {
       $line = $temp[$i];
	$line{$i}=$line;
 	@entry = split(/\t/,$line);
	$queryID_current = $entry[0];
	$bitscore_entry =$entry[11];
	if ($queryID_current ne $queryID_previous) {
		if ($flag==0) {
			$bitscore{$i}=$bitscore_entry;
			$queryID_previous = $queryID_current;
			$flag=1;
		}
		elsif ($flag==1) {
			foreach $key (sort{$bitscore{$b}<=>$bitscore{$a}} keys %bitscore) {
				$value = $bitscore{$key};
				push(@values,$value);
				push(@keys,$key);
			}
			$best_bitscore = $values[0];
			if ($best_bitscore >= 100) {##### modify bitscore cutoff here #####
				$values_index = $keys[0];
				$x = $values_index-$i;
				$line_index = $i + $x;
				$best_line = $line{$line_index};
				print "$best_line\n";
			}
			@values=();
			@keys=();
			%bitscore=();
			$bitscore{$i}=$bitscore_entry;
			$queryID_previous = $queryID_current;
		}
	}
	elsif ($queryID_current eq $queryID_previous) {
		$bitscore{$i}=$bitscore_entry;
		$queryID_previous = $queryID_current;
	}
	$i++;
}

#don't forget to parse last set of entries!
foreach $key (sort{$bitscore{$b}<=>$bitscore{$a}} keys %bitscore) {
	$value = $bitscore{$key};
	push(@values,$value);
	push(@keys,$key);
}
$best_bitscore = $values[0];
if ($best_bitscore >= 100) {##### modify bitscore cutoff here #####
	$values_index = $keys[0];
	$x = $values_index-$i;
	$line_index = $i + $x;
	$best_line = $line{$line_index};
	print "$best_line\n";
}

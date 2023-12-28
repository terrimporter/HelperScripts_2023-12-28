#!/usr/bin/perl
# Teresita Porter, Jan. 21/21
# Script to remove singletons and doubletons from ESV.table as well as Zotus with 0 counts to match results.csv
# USAGE perl filter_ESV_table.plx ESV.table

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $outfile = "ESV.table.parsed";
my $zotu;
my $minval = 3; # minimum cluster size to retain, ex. minval 3 would filter out singletons and doubletons
my $newline;
my $j=0;
my $sum;

# declare arrays
my @in;
my @line;
my @newline;

# declare hashes


# read in ESV.table
open (IN, "<", $ARGV[0]) || die "Error cannot open ESV.table: $!\n";
@in = <IN>; 
close IN;

# create new ESV.table.parsed outfile
open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

# parse ESV.table
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($i==0) { # headers
		print OUT $line."\n";		
	}
	else {
		@line = split(/\t/,$line);	
		$zotu = shift(@line);
		
		# filter my minimum cluster size
		foreach $line (@line) {
			if ($line < $minval) {
				$line = 0;
				$newline[$j] = $line;
				$j++;
			}
			else {
				$newline[$j] = $line;
				$j++;
			}
		}
		$j=0;

		# sum cluster counts
		foreach $newline (@newline) {
			$sum += $newline;
		}
		
		if ($sum == 0) {
			$i++;
			@newline=();
			next;
		}
		else {
			$newline = join "\t", @newline;
			print OUT $zotu."\t".$newline."\n";
		}
	}
	@newline=();
	$sum=0;
	$i++;
}
$i=0;
close OUT;

#!/usr/bin/perl
#April 25, 2014 by Terri Porter
#Script to take OTU size (in header) and expand BLAST results before MEGAN parsing
#Do this so that Nagissa can have stacked bar charts that reflect read abundance instead of just OTU abundance
#USAGE perl expand_blast_results.plx file.blastn

use strict;
use warnings;

#declare var
my $scalar;
my $maxlines;
my $infile;
my $outfile;
my $line;
my $i=0;
my $flag=0;
my $maxsize; #cluster size
my $line2;
my $k=1;
my $scalar2;

#declare array
my @in;
my @record;

open (IN, "<", $ARGV[0]) || die "Error cannot open $ARGV[0]: $!\n";
@in = <IN>;
close IN;

$scalar = scalar(@in);
$maxlines = $scalar-1;

$infile = $ARGV[0];
chomp $infile;

$outfile = $infile.".expanded";

open (OUT, ">>", $outfile ) || die "Error cannot open $outfile:$!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line !~ /Query=/ && $flag == 0) {
		print OUT $line."\n";
	}
	elsif ($line =~ /Query=/ && $flag == 0) {
		$flag = 1;
		$i--;
	}
	elsif ($line =~ /Query=/ && $flag == 1) {
		push(@record, $line);
		if ($line =~ /size_\d+/) {
			$line =~ /size_(\d+)/;
			$maxsize = $1;
		}
		$flag = 2;
	}
	elsif ($line !~ /Query=/ && $flag == 2 ) {
		push(@record, $line);
	}
	elsif ($line =~ /Query=/ && $flag == 2) {
		$flag = 1;

		while ($k <= $maxsize) {
			foreach $line2 (@record) {
				chomp $line2;
				print OUT $line2."\n";
			}
			$k++;
		}
		$k=1;
		@record=();
		$i--;
	}
	$i++;
}
$i=0;

#print last set of records
while ($k <= $maxsize) {
	foreach $line2(@record) {
		chomp $line2;
		print OUT "$line2\n";
	}
	$k++;
}
$k=1;

close OUT;

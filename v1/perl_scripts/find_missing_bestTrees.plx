#!/usr/bin/perl
#March 19, 2013 by Terri Porter
#Script to find missing bestTrees so that the original files so that raxml can be re-run 
#USAGE perl find_missing_bestTrees.plx

use strict;
use warnings;

#declare var
my $total;
my $path;
my $i=0;
my $filename='';
my $output;

print "How many trees in total should there be?:\n";
$total = <STDIN>;
chomp $total;

print "Enter path to bestTrees directory including final backslash:\n";
$path = <STDIN>;
chomp $path;

open (OUT, ">>", "missing.txt") || die "Error could not open missing.txt:$!\n";

while ($i < $total) {
	$filename = "RAxML_bestTree.concatenated19.phy.BS".$i.".out";
	$output = qx( ls $path | grep $filename );

	if (defined $output && length $output == 0) {
		print OUT "$filename\n";
	}
	$i++;
	$output='';
	$filename='';
}
$i=0;
close OUT;

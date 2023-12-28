#!/usr/bin/perl
#June 18, 2012 by Terri Porter
#Script to grep id listed in 12F.contam.seed against 12F.uc (BC1 or BC4) to find 'C' contig line and third field with cluster size
#Usage perl grep_uc.plx 12F.uc
######BE SURE TO MODIFY SEARCH ^D OR ^C IF LOOKING FOR CONTAM OR REGULAR CLUSTERS!!!#####

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $filename;
my $i=0;
my $line;
my $seed;
#my $id;
my $j=0;
my $count;
my $clusterSize;
my $sum;

#declare array
#my @seed;
#my @id;
my @output;
my @line;
my @clusterSize;
my @sorted;
my @sorted2;
my @in;

#declare hash
my %histD;
my %histC;

$filename = $ARGV[0];
chomp $filename;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

while ($in[$j]) {
	$line = $in[$j];
	chomp $line;

	if ($line =~ /^D\t/) {
		#print "$line\n";
		@line = split(/\t/, $line);
		$clusterSize = $line[2];
		$clusterSize = $clusterSize-1;#remove libseed 'L' from clusterSize 'D'
		$histD{$j} = $clusterSize;
		print $clusterSize."\n";
	}
	elsif ($line =~ /^C\t/) {
		#print "$line\n";
		@line = split(/\t/, $line);
		$clusterSize = $line[2];
#		$clusterSize = $clusterSize; #keep seed 'S' in clusterSize
		$histC{$j} = $clusterSize;
		print $clusterSize."\n";
	}
	$j++;
	$line=();
	@line=();
	$clusterSize = ();
}
$j=0;

#print sorted list 
open (OUT, ">>", "histC.txt") || die "Error cannot open outfile1: $!\n";
#@sorted = sort { $histC{$a} cmp $histC{$b} } keys %histC; #sort hash by descending value, get a list of keys
while (($j, $clusterSize) = each (%histC)) {
	print OUT $clusterSize."\n";
}
close OUT;

open (OUT2, ">>", "histD.txt") || die "Error cannot open outfile2: $!\n";
while (($j, $clusterSize) = each (%histD)) {
	print OUT2 $clusterSize."\n";
}
close OUT2;

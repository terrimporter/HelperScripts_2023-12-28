#!/usr/bin/perl
#June 18, 2012 by Terri Porter
#Script to grep id listed in 12F.contam.seed against 12F.uc (BC1 or BC4) to find 'C' contig line and third field with cluster size
#Usage perl grep_uc.plx 12F.uc 12F.contam.seed
######BE SURE TO MODIFY SEARCH ^D OR ^C IF LOOKING FOR CONTAM OR REGULAR CLUSTERS!!!#####

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $filename;
my $i=0;
my $line;
my $seed;
my $id;
my $j=0;
my $count;
my $clusterSize;
my $sum;
my $count2=0;

#declare array
my @seed;
my @id;
my @output;
my @line;
my @clusterSize;

$filename = $ARGV[0];
chomp $filename;

open (SEED, "<", $ARGV[1]) || die "Error cannot open seed file: $!\n";
@seed = <SEED>;
close SEED;

#parse out id's from seed file
while ($seed[$i]) {
	$line = $seed[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w+)/;
		$id = $1;
		#print "id $id\n";
		push(@id, $id);
	}
	$i++;
	$line=();
	$id=();
}
$i=0;

#grep ids in uc and grab cluster size
$count = count (@id); #searching for contam OTUs

while ($id[$i]) {
	$id = $id[$i];
	#print "id $id filename $filename\n";
	@output = qx(grep $id $filename);

	while ($output[$j]) {
		$line = $output[$j];
		chomp $line;

		if ($line =~ /^D/) {
			@line = split(/\t/, $line);
			$clusterSize = $line[2];
			$clusterSize = $clusterSize-1;#remove libseed 'L' from clusterSize 'D'
			push(@clusterSize, $clusterSize);
			$count2++;
		}
		$j++;
		$line=();
		@line=();
		$clusterSize = ();
	}
	$j=0;

	$i++;
	$id=();
	@output=();
}
$i=0;

#sum cluster sizes to figure out total number of contam seqs removed
$sum = sum (@clusterSize);

print "Looking for $count contam OTUs.  Found $count2 contam OTUs containing $sum reads\n";

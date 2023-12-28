#!/usr/bin/perl
#Sept. 6, 2012 by Terri Porter
#Script to parse concatenated raxml_info.tre files
#usage perl parse_raxml_info.plx summary_info_BS0_99.tre

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $treeNum;
my $flag=0;
my $treeScore;

#declare array
my @in;

#declare hash
my %treeNum_treeScore;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>; 
close IN;

#parse the concatenated RAxML_info.tre file
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0 && $line =~ /^raxmlHPC -f d -s/) {
		if ($line =~ /concatenated_final\.phy\.BS\d+/) {
			$line =~ /concatenated_final\.phy\.BS(\d+)/;
			$treeNum = $1;
			#print $treeNum."\n"; #test
			$flag=1;
		}
	}
	elsif ($flag==1 && $line =~ /^Final GAMMA-based Score of best tree -\d+\.\d+/) {
		$line =~ /^Final GAMMA-based Score of best tree (-\d+\.\d+)/;	
		$treeScore = $1;
		#print $treeScore."\n"; #test
		$treeNum_treeScore{$treeNum} = $treeScore;
		$flag=0;
	}
	$i++;

	$line=();
#	$treeNum=();
#	$treeScore=();
}
$i=0;

#print out tree numbers and corresponding tree scores from hash
while( ($treeNum, $treeScore) = each %treeNum_treeScore) {
	print "$treeNum\t$treeScore\n";
}


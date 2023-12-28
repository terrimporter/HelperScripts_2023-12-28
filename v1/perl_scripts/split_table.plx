#!/usr/bin/perl
#Terri Porter, March 30, 2017
#Script to split a global denoised OTU table into individual sample tables
#USAGE perl split_table.plx cat.table

use strict;
use warnings;

#declare var
my $headerline;
my $i=0;
my $header;
my $outfile;
my $colnum;
my $awk;

#declare array
my @in;
my @headers;

#declare hash

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in=<IN>;
close IN;

$headerline = $in[0];
@headers = split(/\t/,$headerline);

while ($headers[$i]) {
	$header = $headers[$i];
	chomp $header;
	$header = $header."_R1"; ##### CUSTOMIZE HERE #####

	if ($i == 0) {
		$i++;
		next;
	}
	else {
		$colnum = $i+1;
		$awk = `awk -v colnum=$colnum 'BEGIN {FS="\t"} {print \$1,\$colnum}' cat.table`;
		$outfile = $header.".table";	
		open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";
		print OUT $awk;
		close OUT;
		$awk=();
		$outfile=();
		$colnum=();

	}
	$i++;
}
$i=0;

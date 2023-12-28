#!/usr/bin/perl
#Terri Porter, Aug. 30, 2010
#Script to combine the full_hit_read_tables into a single table
#usage perl combine_hit_read_tables.plx infile1 infile2 infile3 infile4

use strict;
use warnings;

#declare variables
my $line;
my $i=1; #line counter
my $ref_ID;
my $read_freq;
my $total_read_freq1;
my $k;
my $total_read_freq2;
my $total_read_freq3;
my $total_read_freq4;
my $j=0; #counter

#declare array
my @line;
my @ref_ID1;
my @read_freq1;
my @read_freq2;
my @read_freq3;
my @read_freq4;

open (IN1, '<', $ARGV[0]) || die ("Error: $!\n");

while (<IN1>) {
	$line = $_;
	chomp $line;
	if ($i==1) {
		$i++;
		next;
	}
	elsif ($i > 1) {
		@line = split (/\t/, $line);
		$ref_ID = $line[0];
		push (@ref_ID1, $ref_ID);
		$read_freq = $line[1];
		push (@read_freq1, $read_freq);	
		$i++;
		next;
	}
}
close IN1;
foreach $k (@read_freq1) {
	$total_read_freq1 += $k;
}
#print $total_read_freq1."\n";#test
$i=1;

open (IN2, '<', $ARGV[1]) || die ("Error: $!\n");

while (<IN2>) {
	$line = $_;
	chomp $line;
	if ($i == 1 ) {
		$i++;
		next;
	}
	elsif ($i > 1) {
		@line = split (/\t/, $line);
		$read_freq = $line[1];
		push (@read_freq2, $read_freq);
		$i++;
		next;
	}
}
close IN2;
foreach $k (@read_freq2) {
	$total_read_freq2 += $k;
}
#print $total_read_freq2."\n"; #test
$i=1;

open (IN3, '<', $ARGV[2]) || die ("Error: $!\n");

while (<IN3>) {
	$line = $_;
	chomp $line;
	if ($i==1) {
		$i++;
		next;
	}
	elsif ($i>1) {
		@line = split (/\t/,$line);
		$read_freq = $line[1];
		push (@read_freq3, $read_freq);
		$i++;
		next;
	}
}
close IN3;
foreach $k (@read_freq3) {
	$total_read_freq3 += $k;
}
#print $total_read_freq3."\n"; #test
$i=1;

open (IN4, '<', $ARGV[3]) || die ("Error: $!\n");

while (<IN4>) {
	$line = $_;
	chomp $line;
	if ($i==1){
		$i++;
		next;
	}
	elsif ($i > 1){
		@line = split (/\t/, $line);
		$read_freq = $line[1];
		push (@read_freq4, $read_freq);
		$i++;
		next;
	}
}
close IN4;
foreach $k (@read_freq4) {
	$total_read_freq4 += $k;
}
#print $total_read_freq4."\n"; #test

open (OUT, '>>', "combined_read_freq_table.txt") || die ("Error: $!\n");
print OUT "Ref_ID\tNum_reads_1\tNum_reads_2\tNum_reads_3\tNum_reads_4\n";
while ($ref_ID1[$j]){
	print OUT "$ref_ID1[$j]\t$read_freq1[$j]\t$read_freq2[$j]\t$read_freq3[$j]\t$read_freq4[$j]\n";
	$j++;
}
print OUT "Total_reads\t$total_read_freq1\t$total_read_freq2\t$total_read_freq3\t$total_read_freq4\n";
close OUT;

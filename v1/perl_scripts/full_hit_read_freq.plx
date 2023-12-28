#!/usr/bin/perl
#Terri Porter, Aug.27, 2010
#Script to print list of all referenceIDs and corresponding read frequency if available
#usage $perl full_hit_read_freq.plx reference.fasta hit_read_freq_table.txt

use strict;
use warnings;

#declare variables
my $line;
my $ref_ID;
my $i=0; #index
my $current_ref;
my $read_freq;
my $found_flag=0;

#declare array
my @ref_IDS;
my @line;

open (IN1, '<', $ARGV[0]) || die ("Error: $!\n");

while (<IN1>) {
	$line = $_;
	chomp $line;
	if ($line =~ />\S+/) {
		$line =~ />(\S+)/;
		$ref_ID = $1;
		push (@ref_IDS,$ref_ID);
	}
	else {
		next;
	}
}
close IN1;

open (OUT, '>>', "full_hit_read_freq_table.txt") || die ("Error: $!\n");
print OUT "ReferenceID\tReadFrequency\n";

while ($ref_IDS[$i]){
	$current_ref = $ref_IDS[$i];
	open (IN2, '<', $ARGV[1]) || die ("Error: $!\n");

	while (<IN2>){
		$line = $_;
		chomp $line;
		if ($line =~ /$current_ref/) {
			@line = split (/\t/,$line);
			$read_freq = $line[1];
			$found_flag = 1;
		}
		else {
			next;
		}
	}
	if ($found_flag == 1) {
		print OUT "$current_ref\t$read_freq\n";
	}
	else {
		print OUT "$current_ref\t0\n";
	}
	$i++;
	$found_flag=0;
	close IN2;
}

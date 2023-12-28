#!/usr/bin/perl
#Nov.12,2010 by Terri Porter
#Script to read agrep results from two files and write a non-redundant fasta formatted file
#usage $perl remove_duplicates.pl file.agrep file2.agrep

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $seq;
my $id;
my $j=0;
my $k=0;
my $l=0;
my $m=0;
my $current_id;
my $current_seq;
my $checkid;

#declare arrays
my @file1;
my @file2;
my @line;
my @id;
my @seq;
my @id2;
my @seq2;

open (IN1,"<",$ARGV[0]) || die "Error: $!\n";
@file1 = <IN1>;
close IN1;

open (IN2,"<",$ARGV[1]) || die "Error: $!\n";
@file2 = <IN2>;
close IN2;

while ($file1[$i]) {
	$line = $file1[$i];
	chomp $line;
	@line = split(/\|/,$line);
	$seq = $line[0];
	$id = $line[1];
	push(@id,$id);
	push(@seq,$seq);
	$i++;
}

while ($file2[$j]) {
	$line = $file2[$j];
	chomp $line;
	@line = split(/\|/,$line);
	$seq = $line[0];
	$id = $line[1];
	push(@id2,$id);
	push(@seq2,$seq);
	$j++;
}

open (OUT,">>","agrep.fasta") || die "Error: $!\n";

while ($id[$k]) {
	$current_id = $id[$k];
	$current_seq = $seq[$k];
	print OUT ">$current_id\n$current_seq\n";
	while ($id2[$l]) {
		$checkid = $id2[$l];
		if ($current_id eq $checkid) {
			delete $id2[$l];
			delete $seq2[$l];
		}
		$l++;
	}
	$l=0;
	$k++;
}

while ($id2[$m]) {
	$current_id = $id2[$m];
	$current_seq = $seq2[$m];
	print OUT ">$current_id\n$current_seq\n";
	$m++;
}

#!/usr/bin/perl
#Nov.12,2010 by Terri Porter
#Script to read agrep results from two files and write a non-redundant fasta formatted file
#usage $perl remove_agrep_duplicates.pl file.agrep file2.agrep
#modified to work on a bunch of agrep.out files that were merged using merge.plx
#usage $perl remove_agrep_duplciates_merged.plx file.merged

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $seq;
my $id;
my $j=0;
my $k=0;
my $nonredundant_id;
my $checkid;
my $sequence;
my $flag=0;

#declare arrays
my @file1;
my @line;
my @id;
my @seq;
my @id_nonredundant;

#declare hash
my %ID;

open (IN1,"<",$ARGV[0]) || die "Error: $!\n";
@file1 = <IN1>;
close IN1;

#get all ids and seqs from merged agrep.out files

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

#dereplicate ids and corresponding seqs

%ID = map {$_,1} @id;
@id_nonredundant = keys %ID;

#grab first id & sequence for each non-redundant id

open (OUT,">>","agrep.fasta") || die "Error: $!\n";

while ($id_nonredundant[$j]) {
	$nonredundant_id = $id_nonredundant[$j];

	while ($id[$k]) {
		$checkid = $id[$k];
		$sequence = $seq[$k];
		if ($nonredundant_id eq $checkid) {
			if ($flag==0) {
				print OUT ">$checkid\n$sequence\n";
				$flag=1;
			}
			elsif ($flag==1) {
				$k++;
				next;
			}
		}
		$k++;
	}
	$flag=0;
	$j++;
	$k=0;
}	


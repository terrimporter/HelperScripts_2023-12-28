#!/usr/bin/perl
#Nov.25,2010 by Terri Porter
#Script to check a parsed fasta file with 15-digit alphanumeric ids for duplicates
#usage $perl check_fasta_for_duplicate_ids.plx file.fasta

use strict;
use warnings;

#var
my $line;
my $header;
my $id;
my $i=0;
my $j=0;
my $k=0;
my $l=0;
my $nonredundant_id;

#array
my @ids;
my @ids_nonredundant;
my @file;

#hash
my %ID;

open (IN,"<",$ARGV[0]) || die ("Error: $!\n");
@file = <IN>;
close IN;
#print "@file\n";#test
while ($file[$i]) {
	$line = $file[$i];
	chomp $line;
	#print "$line\n";#test
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /(\w{14})/;
		$id = $1;
		push (@ids, $id);
	}
	$i++;
}

#dereplicate

%ID = map {$_,1} @ids;
@ids_nonredundant = keys %ID;

open (OUT,">>","redundancy_check.txt") || die ("Error: $!\n");

while ($ids_nonredundant[$j]) {
	$nonredundant_id = $ids_nonredundant[$j];
	
	while ($ids[$k]) {
		$id = $ids[$k];
		if ($nonredundant_id eq $id) {
			$l++;
		}
		$k++;
	}
	print OUT "$nonredundant_id\t$l\n";
	$l=0;
	$k=0;
	$j++;
}
		

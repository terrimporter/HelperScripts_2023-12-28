#!/usr/bin/perl
#Terri Porter, Aug.26, 2010
#Script to filter UBLAST output according to the highest %identity per query sequence
#usage $perl parse_ublast.plx infile.b6

use strict;
use warnings;

#declare variables
my $line;
my $query_long;
my $query;
my $c; #element
my $f; #element
my $bitscore;
my $m=0; #indexing
my $best_bitscore;
my $i=0;
my $q;
my $bi;
my $flag=0;
my $k=0; #indexing
my $j=0; #counter
my $first_field;
my $query_ID;
my $hit_ID;
my $l=0; #indexing

#declare arrays
my @line;
my @query_IDS;
my @query_IDS_nonredundant;
my @file;
my @bitscore;
my @sorted_bitscore;
my @best_bitscore;
my @file2;

#declare hash
my %query_IDS;

get_all_query_IDs();
dereplicate_array_of_ids();
find_best_bitscore();
print_top_hit_table();
print_query_hit_map();

#subroutine to get all query IDs

sub get_all_query_IDs {

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /\d+|\S+|\w{14}\s+/) {
		@line = split(/\t/,$line);
		$query_long = $line[0];
		if ($query_long =~ /\d+|\S+|\w{14}\s+/) {
			$query_long =~ /(\w{14})\s+/;
			$query = $1;
			push (@query_IDS, $query);
		}
	}
}
close IN;
print "got all query IDS\n"; #status message
}

#subroutine to dereplicate array of ids

sub dereplicate_array_of_ids {

%query_IDS = map {$_, 1} @query_IDS;
@query_IDS_nonredundant = keys %query_IDS;

print "dereplicated array of query IDs\n"; #status message
}

#subroutine to get best bitscore
sub find_best_bitscore {

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
@file = <IN>;
close IN;
print "read file into array and searching for best bitscore\n"; #status

while ($query_IDS_nonredundant[$i]) {
	$c = $query_IDS_nonredundant[$i];
	#print "$c\n";#test

	while ($file[$m]) {
		#print "started reading file in array\n"; #test
		$f = $file[$m];
		if ($f =~ /$c/) {
			#print "\t\tfinding matches\n"; #test
			@line = split(/\t/, $f);
			$bitscore = $line[11];
			push (@bitscore,$bitscore);
		}
		else {
			#print "\tnot finding matches\n"; #test
		}
		$m++;
		#print "."; #status message
	}
	#print "\n";

	@sorted_bitscore = sort { $a<=> $b } @bitscore;
	#print "@sorted_bitscore\n"; #test
	$best_bitscore = pop(@sorted_bitscore);
	#print "$best_bitscore\n"; #test
	push (@best_bitscore, $best_bitscore);
	@bitscore=();
	@sorted_bitscore=();
	$i++;
	$m=0;
}
print "got best bitscores\n"; #status
}

#subroutine to find line to print top hit table

sub print_top_hit_table {
open (OUT1, '>>', "top_hit.table") || die ("Error: $!\n");

while ($query_IDS_nonredundant[$j]) {
	$q = $query_IDS_nonredundant[$j];
	$bi = $best_bitscore[$j];
	while ($file[$k]) {
		$f = $file[$k];
		if ($flag==0 && $f =~ /$q.+$bi/) {
			print OUT1 $f."\n";
			$flag=1; #if two lines with same query and bitscore, print the first one only
		}
		$k++;
	}
	$flag=0;
	$k=0;
	$j++;
}
close OUT1;
print "printed top hit table\n"; #status
}

#subroutine to print list of query IDs and top hit reference IDs to be used for printing read frequency

sub print_query_hit_map {

open (IN2, '<', "top_hit.table") || die ("Error: $!\n");

open (OUT2, '>>', "query_hit.table" ) || die ("Error: $!\n");
print OUT2 "QueryID\tHitID\n";

while (<IN2>) {
	$line = $_;
	chomp $line;
	if ($line =~ /\w{14}\s{1}.+\tNEWLEP\S+\t/) {
		$line =~ /(\w{14}).+\t(NEWLEP\S+)\t/;
		$query_ID = $1;
		$hit_ID = $2;
		print OUT2 "$query_ID\t$hit_ID\n";
	}
}
close IN2;
close OUT2;
print "printed query hit map\n"; #status
}

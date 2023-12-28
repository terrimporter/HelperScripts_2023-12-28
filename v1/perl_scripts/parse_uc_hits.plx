#!/usr/bin/perl
#Terri Porter, Aug.26, 2010
#Script to parse results of query search aganist reference library: hits (H), no-hits (N)
#usage parse_uc_hits.plx infile.uc

use strict;
use warnings;

#declare variables
my $line;
my $query;
my $ref;
my $cluster_number;
my $i=0; #counter
my $to_match;
my $j=0; #counter

#declare arrays
my @line;
my @clusters;
my @refs;

print_query_results();
parse_ref_ids_and_cluster_numbers();
count_number_hits_per_cluster();

#subroutine to summarize query results for hits and no-hits

sub print_query_results {

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
open (OUT1, '>>', "hit_summary.txt") || die ("Error: $!\n");
open (OUT2, '>>', "nohit_summary.txt") || die ("Error: $!\n");

print OUT1 "query\tref\n";
print OUT2 "query\tref\n";

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^H/) {
		@line = split (/\t/,$line);
		$query = $line[8];
		$ref = $line[9];
		print OUT1 "$query\t$ref\n";
	}
	elsif ($line =~ /^N/) {
		@line = split (/\t/,$line);
		$query = $line[8];
		$ref = $line[9];
		print OUT2 "$query\t$ref\n";
	}
}
close IN;
close OUT1;
close OUT2;
}

#subroutine to get referenceID and corresponding cluster numbers

sub parse_ref_ids_and_cluster_numbers {

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^L/) {
		@line = split (/\t/,$line);
		$cluster_number = $line[1];
		push (@clusters, $cluster_number);
		$ref = $line[8];
		push (@refs, $ref);
		#print OUT3 "$ref\n";
	}
}
#print "@clusters\t"; #test
#print "@refs\t"; #test
close IN;
#close OUT3;
}

#subroutine to count number of hits per cluster

sub count_number_hits_per_cluster {

open (OUT3, '>>', "refs_hit_summary.txt") || die ("Error: $!\n");
print OUT3 "cluster_numer\treference_ID\tnumber_hits\n";

while ($clusters[$i]) {
	$to_match = $clusters[$i];
	#print OUT3 "$to_match\ttomatchvar\t"; #test
	open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^H/) {
			@line = split (/\t/, $line);
			$cluster_number = $line[1];
			#		print OUT3 "$cluster_number\tclusternumbervar\n"; #test
			if ($to_match == $cluster_number) {
				$j++;
			}
			else {
				next;
			}
		}
	}
	close IN;
	$ref = $refs[$i];
	print OUT3 "$to_match\t$ref\t$j\n";
	$j=0;
	$i++;	
}
close OUT3;

}



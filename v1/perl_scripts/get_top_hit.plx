#!/usr/bin/perl
#Terri Porter, Aug.31, 2010
#Script to parse through blastall tabular output format -m 8 (no comment lines) and grab just the top hit for each unique queryID
#usage perl get_top_hit.plx in

use strict;
use warnings;

#declare variables
my $line;
my $queryID;
my $i=0; #indexing
my $j=1; #counter
my $query_ref;
my $query_to_match;
my $entry;

#declare arrays
my @line;
my @queryIDs;
my @query_IDs_nonredundant;
my @keep;

#declare hash 
my %queryIDs;

get_query_IDs();
dereplicate_query_IDs();
get_top_hit();

#subroutine to grab query IDS from blast output

sub get_query_IDs {

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
	@line = split (/\t/,$line);
	$queryID = $line[0];
	push (@queryIDs, $queryID);
}
close IN;

}
#print "@queryIDs\n"; #test

#subroutine to dereplicate array of query IDs

sub dereplicate_query_IDs {

%queryIDs = map {$_, 1} @queryIDs;
@query_IDs_nonredundant = keys %queryIDs;
#print "@query_IDs_nonredundant\n";
}

#subroutine to grab the first hit for each queryID

sub get_top_hit {

open (OUT, '>>', "top_blast_hits.txt") || die ("Error: $!\n");

while ($query_IDs_nonredundant[$i]) {
	$query_ref = $query_IDs_nonredundant[$i];
	open (IN, '<', $ARGV[0]) || die ("Error: $!\n");

	while (<IN>) {
		$line = $_;
		chomp $line;
		$entry = $line;
		@line = split (/\t/,$line);
		$query_to_match = $line[0];
		if ($query_ref eq $query_to_match) {
			push (@keep, $entry);
			#print "@keep\n"; #test
		}
		else {
			next;
		}
	}
	print OUT $keep[0]."\n";
	@keep=();
	$i++;
	close IN;
}
close OUT;
}

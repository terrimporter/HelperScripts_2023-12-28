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
my $i=0; #counter
my $current_query;
my $bitscore;
my $best_bitscore;
my $hit;
my $j=0; #counter

#declare arrays
my @line;
my @query_IDS;
my @query_IDS_nonredundant;
my @bitscore;
my @hits;
my @sorted_bitscore;
my @best_bitscores;

#declare hash
my %query_IDS;

get_all_query_IDs();
dereplicate_array_of_ids();
find_best_bitscore();
print_best_hit_summary();

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
	else {
		next;
	}
}
close IN;
}

#subroutine to dereplicate array of ids

sub dereplicate_array_of_ids {

%query_IDS = map {$_, 1} @query_IDS;
@query_IDS_nonredundant = keys %query_IDS;

}

#subroutine to get best bitscore
sub find_best_bitscore {
while ($query_IDS_nonredundant[$i]) {
	$current_query = $query_IDS_nonredundant[$i];
	open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /$current_query/) {
			@line = split(/\t/, $line);
			$hit = $line[1];
			$bitscore = $line[11];
			push (@bitscore,$bitscore);
			push (@hits,$hit);
		}
		else {
			next;
		}
	}
	close IN;
	@sorted_bitscore = sort { $a<=> $b } @bitscore;
	$best_bitscore = pop(@sorted_bitscore);
	push (@best_bitscores, $best_bitscore);
	@bitscore = ();
	@sorted_bitscore = ();
	$i++;
}
}

#subroutine to print summary of top hits

sub print_best_hit_summary {

open (OUT,'>>',"best_hit_summary") || die ("Error: $!\n");

print OUT "QueryID\tHitID\tBitscore\n";
while ($query_IDS_nonredundant[$j]) {
	print OUT "$query_IDS_nonredundant[$j]\t$hits[$j]\t$best_bitscores[$j]\n";
	$j++;
}
close OUT;
}

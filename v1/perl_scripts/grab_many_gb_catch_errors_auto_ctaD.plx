#!/usr/bin/perl
#April 12, 2018 add ctaD[GENE] to entrez query to find more bacterial cytochrome c oxidase sequences
#March 8, 2013 edited to work with lists of 100 taxa names at a time, and auto rerun
#March 7, 2013 by Terri Porter
#script to replace use of ebot_nucleotide.plx and filter_results_numberingfixed.plx and filter_results_pickup.plx
#will need to add another script to do something similar to quick_filter_gb_results.plx
#Written to catch errors and move on instead of crashing the whole f-ing job
#USAGE perl grab_many_gb_catch_errors.plx 100_foramtted_species.txt

use strict;
use warnings;
use Bio::DB::EUtilities;

my $factory;
my $count;
my $hist;
my $retry;
my $out;
my $taxonlist;
my $filename;
my $term;
my @in;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
$taxonlist = $in[0];
chomp $taxonlist;

$term = "\"ctaD\"[gene] AND (".$taxonlist.")";

$factory = Bio::DB::EUtilities -> new (-eutil => 'esearch',
										-email => 'terriblue2002@yahoo.com',
										-db => 'nucleotide',
										-term => $term,
										-usehistory => 'y');

$count = $factory -> get_count;

#get history from queue
$hist = $factory -> next_History || die 'No history data returned';
print "History returned\n";

#db carries over from above
$factory -> set_parameters (-eutil => 'efetch',
							-rettype => 'gb',
							-history => $hist);

$retry = 0;
my ($retmax, $retstart) = (500,0);

$filename = $ARGV[0];
$filename =~ s/\.txt//;
$filename = $filename."_seqs.gb";

open ($out, ">", $filename) || die "Can't open file: $!\n";

RETRIEVE_SEQS:
while ($retstart < $count) {
	$factory -> set_parameters (-retmax => $retmax,
								-retstart => $retstart);
	
	eval {
		$factory -> get_Response (-cb => sub {
											my ($data) = @_;
											print $out $data
										 }
								 );
	};

	if ($@) {
		die "Server error : $@.  Try again later" if $retry == 5;
		print STDERR "Server error, redo #$retry\n";
		$retry++ && redo RETRIEVE_SEQS;
	}

	print "Retrieved $retstart\n";
	$retstart += $retmax;

}

close $out;


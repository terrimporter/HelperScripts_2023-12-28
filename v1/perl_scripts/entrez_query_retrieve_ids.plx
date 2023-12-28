#!/usr/bin/perl
# Teresita M. Porter, July 3/20
# Script to submit an Entrez query to nucleotide database, retrieve ids for later retrieval
# narrow search by SLEN and PDAT (publication date) to reduce records returned

# https://bioperl.org/howtos/EUtilities_Cookbook_HOWTO.html#item23

use Bio::DB::EUtilities;
use strict;
use warnings;

# vars
my $query;
my $out;
my $filename;
my $i=0;

# array
 my @years = ("1982:1998","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020");
# my @years = ("2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020");
# my @years = ("1982:1998");

# loop through years and do entrez search
while ($years[$i]) {
	$query = $years[$i];
	chomp $query;
	$filename = $query.".gb";

	submit_entrez_query($query,$filename);

	$i++;
}
$i=0;


########


sub submit_entrez_query{

	$query=$_[0];
	$filename=$_[1];

	my $term;
	my $factory;

	$term = '(COI[gene] OR CO1[gene] OR coxI[gene] OR cox1[gene]) AND 500:1500[SLEN] AND ('.$query.'[PDAT])';

	# set optional history queue
	$factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
                                       -email      => 'terriblue2002@yahoo.com',
                                       -db         => 'nucleotide',
                                       -term       => $term,
                                       -usehistory => 'y');
	my $count = $factory->get_count;

	# get history from queue
	my $hist = $factory->next_History || die 'No history data returned';
	print "History returned\n";

	# note db carries over from above
	# use rettype 'gbwithparts' instead of 'gb' so that CONTIG records are returned in full
	$factory->set_parameters(-eutil => 'efetch',
                         -rettype => 'gbwithparts',
                         -history => $hist);

	my $retry = 0; 
	my ($retmax, $retstart) = (500,0);

	open ($out, ">", $filename) || die "Canot open outfile: $!\n";

	# retrieve seqs subroutine that is called implicitly
	RETRIEVE_SEQS:
	while ($retstart < $count) {
    	$factory->set_parameters(-retmax => $retmax,
        	                     -retstart => $retstart);

	    eval{
    	    $factory->get_Response(-cb =>
        	    sub {my ($data) = @_; print $out $data} );
	    };
    	if ($@) {
        	die "Server error: $@.  Try again later" if $retry == 5;
 	       print STDERR "Server error, redo #$retry\n";
    	    $retry++ && redo RETRIEVE_SEQS;
	    }
    	print "Retrieved $retstart";
    	$retstart += $retmax;
	}

	close $out;

}

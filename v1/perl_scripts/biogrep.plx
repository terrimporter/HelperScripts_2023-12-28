#!/usr/bin/perl
#Nov.1,2010 by Terri Porter
#Script to use Bio::Grep module to create an agrep database, search database using mid primer sequences, use this output with biogrep_primersearch.plx
#usage $perl biogrep_midsearh.plx mid.txt

use strict;
use warnings;
use Bio::Grep;
use Bio::SeqIO;

#declare variables
my $db_object;
my $search_result;
my $primer;
my $id;
my $query;
my $i=0;
my $filename;
my $string;
my $string2;
my $found_id;

#declare arrays
my @found_ids;
my @id;
my @primer;
my @primer2;
my @query;

#open (OUT,">>", "biogrep.out") || die ("Error: $!\n");

open (PRIMER, "<", $ARGV[0]) || die ("Error:$!\n");
@primer = <PRIMER>;
close PRIMER;

grab_query_strings();
create_agrep_dbase();
while ($query[$i]) {
	search_agrep_dbase();
	$i++;
}

#grab primer sequences to be used as query strings for agrep search
sub grab_query_strings {

foreach $primer (@primer) {
	chomp $primer;
	@primer2 = split (/\t/, $primer);
	$id = $primer2[0];
	push (@id,$id);
	$query = $primer2[1];
	push (@query,$query);
}
#print "@query\n"; #test
}

#create a database for Agrep to search, only need to do this once

sub create_agrep_dbase {

$db_object = Bio::Grep->new('Agrep');

$db_object->generate_database( {
		file => '9.TCA.454Reads.fna',
		datapath => '/home/terri/Saina/',
	});

}

#search the Agrep database created above

sub search_agrep_dbase {
$string = $query[$i];
$db_object->search ({
	query => $string,
	mismatches => 1,
	reverse_complement => 0,
	database => '9.TCA.454Reads.fna',
});

#push list of ids into array, then print fasta file for these
while ($search_result = $db_object->next_res ) {
	$filename = $id[$i];
	#open (OUT,">>","$filename.out") || die ("Error:$!\n");
	#print OUT $search_result->sequence->id . "\n";
	push @found_ids, $search_result->sequence_id;
}

my $seq_in = $db_object->get_sequences(\@found_ids);

my $seq_out = Bio::SeqIO->new (-file => ">$filename.out",
	                       -format => 'fasta');

while (my $seq = $seq_in ->next_seq) {
	$seq_out->write_seq($seq);
}
@found_ids=();

}



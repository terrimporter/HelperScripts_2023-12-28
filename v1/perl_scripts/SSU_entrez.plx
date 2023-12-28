#!/usr/bin/perl

#Script to fetch sequences from GenBank using an Entrez query
#September 22, 2009 written by Terri Porter
#Usage $perl fetchseqs.plx > out.txt
#modified to grab SSU and automatically write outfile
#modified usage $perl SSU_entrez.plx
#regular Bio::DB::GenBank can't handle over 37000 seqs, so use EUtilities intead

use warnings;
use strict;

use Bio::DB::EUtilities;
use Perl6::Say;

#declare var
my $factory;
my $count;
my $hist;
my $retry;
my $retmax;
my $retstart;
my $out;
my $data;

$factory = Bio::DB::EUtilities->new(-eutil	=> 'esearch',
					-db	=> 'nucleotide',
					-term   =>'"18S" AND "rRNA"',           
				-usehistory	=> 'y');

$count = $factory->get_count;
$hist = $factory->next_History || die 'No history data returned';
print "History returned\n";

$factory-> set_parameters(	-eutil	=> 'efetch',
				-rettype	=> 'fasta',
				-history	=> $hist);

$retry = 0;

($retmax, $retstart) = (500,0);

open ($out, '>', 'seqs.fa') || die "Can't open file: $!\n";

RETRIEVE_SEQS:

while ($retstart< $count) {
	$factory -> set_parameters(	-retmax	=> $retmax,
					-retstart	=> $retstart);
	eval {
		$factory->get_Response (-cb	=> sub {($data) = @_; print $out $data} );
#		$factory->get_Response(-file	=> sub {$data = @_; print $out $data} );
	};
	if ($@) {
		die "Server error: $@.  Try again later" if $retry == 5;
		print "Server error, redo #$retry\n";
		$retry++ && redo RETRIEVE_SEQS;
	}
	say "Retrieved $retstart";
	$retstart += $retmax;
}

close $out;


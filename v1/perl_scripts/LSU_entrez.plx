#!/usr/bin/perl
#Sept. 15, 2011 edited to search for all fungal LSU genbank seqs 
#Script to fetch sequences from GenBank using an Entrez query
#September 22, 2009 written by Terri Porter
#Usage $perl fetchseqs.plx

use warnings;
use strict;

use Bio::DB::GenBank;     
use Bio::DB::Query::GenBank;
use Bio::SeqIO; 
use Bio::DB::EUtilities;
use Statistics::Lite qw(:all);
use Perl6::Say;

#perform entrez query, retrieve fasta
print "\nDoing entrez query...\n";

#var
my $factory;
my $count;
my $hist;
my $retry;
my $retmax;
my $retstart;
my $out;
my $data;

#Bio::DB::GenBank can't handle over 37,000 sequences, so use EUtilities instead see SSU_entrez.plx for an example of how this works

$factory = Bio::DB::EUtilities -> new	(-eutil	=> 'esearch',
					-db	=> 'nucleotide',
					-term => 'Fungi AND ("large ribosomal subunit" OR 28S OR 26S OR 25S) NOT (mitochondrial OR mitochondrion)',
			     		-usehistory	=> 'y');		
$count = $factory -> get_count;
$hist = $factory -> next_History || die 'No history data returned';
print "History returned\n";
$factory -> set_parameters (	-eutil	=> 'efetch',
				-rettype	=> 'fasta',
				-history	=> $hist);
$retry = 0;
($retmax, $retstart) = (500,0);

open ($out, ">", "entrez.fasta") || die "Can't open file: $!\n";

RETRIEVE_SEQS:

while ($retstart< $count) {
	$factory -> set_parameters(	-retmax	=> $retmax,
					-retstart	=> $retstart);
	eval {
		$factory->get_Response (-cb	=> sub {($data) = @_; print $out $data} );
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

print "\nParsing fasta to get gi list...\n";

#var
my $i=0;
my $line;
my $gi;
my $gb;
my $partial;
my $genus;
my $species;

#array
my @in;
my @gi;
my @gb;
my @partial;
my @line;

open (IN,"<","entrez.fasta") || die ("Error cannot read from entrez.fasta: $!\n");
@in = <IN>;
close IN;

open (OUT2, ">>", "gi.query") || die ("Error cannot write to gi.query: $!\n");
open (OUT2b,">>","gb.query") || die ("Error :$!\n");

while($in[$i]) {
	$line = $in[$i];
	chomp $line;
	if ($line =~ /^>/) {
		@line = split (/\|/,$line);
		$gi = $line[1];
		push(@gi,$gi);
		$gb = $line[3];
		push(@gb,$gb);
		#$partial = $line[4];
		#@partial = split(/ /,$partial);
		#$genus = $partial[0];
		#$species = $partial[1];
		print OUT2 "$gi\n";
		print OUT2b "$gb\n";
		#print OUT3 "$gb\n";
	}
	$i++;
}
close OUT2;
#close OUT3;

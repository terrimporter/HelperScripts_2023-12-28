#!/usr/bin/perl

#Script to fetch sequences from GenBank using an Entrez query
#September 13, 2010 written by Terri Porter
#Usage $perl fetchseqs.plx > out.txt

use warnings;
use strict;

use Bio::DB::GenBank;     
use Bio::DB::Query::GenBank;
use Bio::SeqIO; 

my $query = Bio::DB::Query::GenBank->new         
	(-query   	=>'Saccharomyces[ORGN] AND "internal transcribed"',           
	-db      	=> 'nucleotide');

#get a Genbank database handle     
my $gb = Bio::DB::GenBank->new
	(-format	=> 'Fasta');
			 
my $stream = $gb->get_Stream_by_query($query);     

while( my $seqobj =  $stream->next_seq() ) {       
 	my $seq = $seqobj->seq();
	my $ID = $seqobj->id();
	my $desc = $seqobj->desc();
	print ">",$ID,"          ",$desc,"\n",$seq,"\n";	
}     



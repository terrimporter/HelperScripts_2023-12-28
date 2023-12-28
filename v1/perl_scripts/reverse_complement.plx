#!/usr/bin/perl
#Nov. 24, 2016 by Terri Porter
#Script to reverse complement wierd ITS sequences where R1 reads are off the ITS4_F primer and R2 reads are off the ITS1F_R primer.  In this case, need to reverse complement all R1 reads before using usearch uchime2_ref method because it requires plus stranded fastas.
#USAGE perl reverse_complement.plx infile.fasta

use strict;
use warnings;

use Bio::SeqIO;

#declare var
my $seq;
my $revcom;
my $basename;
my $outfilename;
my $in;
my $out;

#declare array
my @filenameparts;

$in = Bio::SeqIO->new (-file => "$ARGV[0]",
						-format => 'fasta');

@filenameparts = split(/\./,$ARGV[0]);
$basename = $filenameparts[0];
$outfilename = $basename.".revcom.fasta";

$out = Bio::SeqIO->new(-file => ">>$outfilename",
					-format => 'fasta');

while ( $seq = $in->next_seq()) {
	$revcom = $seq->revcom;
	$out->write_seq($revcom);
}

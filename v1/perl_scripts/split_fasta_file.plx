#!/usr/bin/perl
#October 9, 2009 by Terri Porter
#Script to split a large fasta file into smaller fasta files.  Necessary to get Chimera Check server to run
#USAGE $perl splitfastafile.plx < bigfasta.txt

use strict;
use warnings;

#declare variables
my $line;
my $num_elements;
my $num_parts;
my $num_parts_roundup;
my $i=1;
my $header_part;
my $seq_part;
my $j=0;

#declare arrays
my @header;
my @seq;
my @header_part;
my @seq_part;

while (<>){
	$line =$_;
	chomp($line);
	if (/>/){
		push (@header, $line);
	}
	else {
		push (@seq, $line);
	}
}

$num_elements = scalar(@header);
$num_parts = $num_elements/50000;#edit size here too
$num_parts_roundup = $num_parts + 1;

while ($i <= $num_parts_roundup){
	@header_part = splice(@header,0,50000);#edit size of files manually if needed
	@seq_part = splice(@seq,0,50000);#edit here too
	
	open (OUT, ">>", "part$i\.out") || die "Can't open 'part.$i\.out': $!";
	
	while ($header_part[$j]){
		print OUT $header_part[$j],"\n";
		print OUT $seq_part[$j],"\n";
		$j++;
	}
	$i++;
	$j=0;#reset counter
	@header_part=();#empty array
	$seq_part = ();
	
	close OUT;
}	


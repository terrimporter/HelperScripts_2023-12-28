#!/usr/bin/perl
#January 27, 2010 by Terri Porter
#Script to split a large ids.infile into smaller files that GenBank will accept when fetching genbank records. 
#USAGE $perl split_ids_outfile.plx < big_ids.outfile

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
my @ids;
my @ids_part;

while (<>){
	$line =$_;
	chomp($line);
	push (@ids, $line);
}

$num_elements = scalar(@ids);
$num_parts = $num_elements/500;#edit size here too
$num_parts_roundup = $num_parts + 1;

while ($i <= $num_parts_roundup){
	@ids_part = splice(@ids,0,500);#edit size of files manually if needed					
	open (OUT, ">>", "part$i\.out") || die "Can't open 'part.$i\.out': $!";
							
	while ($ids_part[$j]){
		print OUT $ids_part[$j],"\n";
		$j++;
	}
	$i++;
	$j=0;#reset counter
	@ids_part=();#empty array								
	close OUT;
}

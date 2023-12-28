#!/usr/bin/perl
#edited Sept.11, 2012 to work with agrep to allow one mismatch for Tara.  
#See part1_primer_sort_010312.plx for how to reformat and use agrep
#Sept.10, 2012 by Terri Porter
#Script to search for primer binding sites in a list of fasta seqs (for Tara)
#primer.txt file should only contain one specific primer or all specific primers for a single degenerate primer sequence
#usage perl primer_search_agrep.plx reference.fasta primer.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $primerName;
my $primerSeq;
my $j=0;
my $line2;
my $accession;
my $k=0;
my $count=0;
my $pattern;

#declare array
my @ref;
my @primer;
my @line;
my @output;

#declare hash
my %primerName_Seq;

open (REF, "<", $ARGV[0]) || die "Error cannot open reference.fasta: $!\n";
@ref = <REF>;
close REF;

open (PRIMER, "<", $ARGV[1]) || die "Error cannot open primer.txt: $!\n";
@primer = <PRIMER>;
close PRIMER;

#parse primer file
while ($primer[$i]) {
	$line = $primer[$i];
	chomp $line;
	
	@line = split(/\t/,$line);
	$primerName = $line[0];
	$primerSeq = $line[1];
	$primerName_Seq{$primerName} = $primerSeq;

	$i++;
	$line=();
	@line=();
	$primerName=();
	$primerSeq=();
}
$i=0;

open (TEMP, ">>", "temp.agrep") || die "Error cannot open temp.agrep: $!\n";

#use grep to search reference fasta file
while ($ref[$i]) {
	$line = $ref[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$accession = $1;
		$j=$i+1;
		$line2 = $ref[$j];
		chomp $line2;
		#print $line2."\n";#test
		
#		while(($primerName,$primerSeq) = each %primerName_Seq) {
			#search for exact primer matches only
#			if ($line2 =~ /$primerSeq/) {
#				$count++;
#				print "found a primer match\n";#test
#			}
#		}

		#reformat for an agrep search
		print TEMP "$line2|$accession\n";
	}
	$i++;
	$line=();
	$accession=();
	$j=0;
	$line2=();
	#print "$count primer matches found\n";
}
$i=0;
close TEMP;

while (($primerName, $primerSeq) = each %primerName_Seq) {
	$pattern = $primerSeq;
#	print "pattern: $pattern\n";#test

	###adjust agrep settings here###
	@output = qx(agrep -i -1 -D1 -I1 "$pattern" temp.agrep); ###change -# to control number of allowed mismatches###
	
	while ($output[$j]) {
		$line = $output[$j];
		print "line: $line\n";#test
		$count++;
		$j++;
	}
	$j=0;
	$line=();
	$pattern=();
	@output=();
}

unlink("temp.agrep");

print "$count degenerate primer matches found\n";

	

#!/usr/bin/perl
#Dec.9, 2016 by Terri Porter
#Script to make a spreadsheet from RDP classifier CO1v2 output
#USAGE perl make_matrices_from_RDP2.plx 2014_mapping.txt YEAR_rdp.out

use strict;
use warnings;

#declare var
my $line;
my $seqID;
my $fieldID;
my $i=0;
my $outfile;
my $header;
my $order;
my $orderBP;
my $family;
my $familyBP;
my $genus;
my $genusBP;

#declare array
my @mapping;
my @line;
my @rdp;
my @header;

#declare hash
my %map; #key=seqID, value=fieldID

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@mapping = <IN>;
close IN;

foreach $line (@mapping){
	chomp $line;
#	print "line: $line\n"; #test
	@line = split(/\t/,$line);
	$seqID = $line[0];
#	print "seqID: $seqID\n"; #test
	$fieldID = $line[1];
#	print "fieldID: $fieldID\n"; #test
	$map{$seqID} = $fieldID;
}
$line=();

open (IN2, "<", $ARGV[1]) || die "Error cannot open infile2: $!\n";
@rdp = <IN2>;
close IN2;

$outfile = $ARGV[1].".txt";

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";
print OUT "SampleID\tFieldSite\tOrder\tOrderBP\tFamily\tFamilyBP\tGenus\tGenusBP\n";

while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;
	print "line: $line\n"; #test

	@line = split(/\t/,$line);
	$header = $line[0];
	@header = split(/;/,$header);
	$seqID = $header[1];
	print "seqID: $seqID\n"; #test
	$order = $line[10];
	$orderBP = $line[11];
	$family = $line[12];
	$familyBP = $line[13];
	$genus = $line[14];
	$genusBP = $line[15];

	if (exists $map{$seqID}) {
		$fieldID = $map{$seqID};
		print OUT "$seqID\t$fieldID\t$order\t$orderBP\t$family\t$familyBP\t$genus\t$genusBP\n";
	}
	else {
		print "Cannot map seqID $seqID to a fieldID because it's not one of the analyzed samples here.\n";
	}

	$i++;
}
$i=0;
close OUT;

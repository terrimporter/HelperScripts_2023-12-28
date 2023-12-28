#!/usr/bin/perl
#Feb.6, 2012 by Terri Porter
#Script to grab fastas from .fasta and .qual files according to a list of read IDs
#usage perl get_fastas_from_list.plx file.fasta OR file.qual AND readID.list (one id per line)

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $header;
my $j;
my $headerID;
my $seq;
my $readID;

#declare array
my @fasta;
my @list;

#declare hash
my %fasta;
my %filteredFasta;
	
open (IN,"<", $ARGV[0]) || die "Error cannot read from infile: $!\n";
@fasta = <IN>;
close IN;

open (IN2,"<",$ARGV[1]) || die "Error cannot read id list: $!\n";
@list = <IN2>;
close IN2;

#put fasta file into hash indexed by readID

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w+)/;
		$headerID = $1;
		$j = $i+1;
		$seq = $fasta[$j];
		chomp $seq;
		$fasta{$headerID} = $seq;
		
		$header=();
		$headerID=();
		$j=();
		$seq=();

	}
	$i++;
}
$i=0;

#grab just the fastas with readIDs of interest

while ($list[$i]) {
	$readID = $list[$i];
	chomp $readID;

	if ($fasta{$readID} ) {
		$seq = $fasta{$readID};
		$filteredFasta{$readID} = $seq;
	}
	$i++;

	$readID=();
	$seq=();

}
$i=0;

#print filtered fastas to a file

open (OUT,">>","outfile.fasta") || die "Error cannot write to outfile.fasta: $!\n";

while (my ($key,$value) = each (%filteredFasta) ) {
	print OUT ">".$key."\n".$value."\n";
}
close OUT;

#!/usr/bin/perl
#Script to build a qual file for marker.nosingletons.fasta after filtering out contaminants
#usage perl get_qual_for_fasta.plx marker.nosingletons.fasta marker.qual.sorted

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $header;
my $id;
my $j;
my $line2;
my $qualSeq;

#declare array
my @fasta;
my @qual;

#declare hash
my %qual;

open (FASTA, "<", $ARGV[0]) || die "Error cannot open fasta file: $!\n";
@fasta = <FASTA>;
close FASTA;

open (QUAL, "<", $ARGV[1]) || die "Error cannot open qual file: $!\n";
@qual = <QUAL>;
close QUAL;

#put qual file into hash
while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/^>//;
		$id = $header;
		$j=$i+1;
		$line2 = $qual[$j];
		chomp $line2;
		$qualSeq = $line2;
		$qual{$id} = $qualSeq;
	}
	$i+=2;
	$line=();
	$header=();
	$id=();
	$j=();
	$line2=();
	$qualSeq=();
}
$i=0;

#parse through fasta file and build new qual file

open (OUT, ">>", "nosingletons.qual") || die "Error cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/^>//;
		$id = $header;

		if (exists $qual{$id}) {
			$qualSeq = $qual{$id};
			print OUT ">$id\n$qualSeq\n";
		}
		else {
			print OUT "ERROR $id\n";
		}
	}
	$i+=2;
	$line=();
	$header=();
	$id=();
	$qualSeq=();
}
$i=0;

close OUT;

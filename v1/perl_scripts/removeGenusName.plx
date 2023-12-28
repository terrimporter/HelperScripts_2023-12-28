#!/usr/bin/perl
#March 21, 2012 by Terri Porter
#Script to remove genus names from FungiLSU_train_1400bp_8506_mod.fasta / mytrainseq.fasta to avoid parsing errors once trained rank set to family
#usage perl removeGenusName.plx FungiLSU_train_taxid.txt/mytaxon.txt FungiLSU_train_1400bp_8506_mod.fasta/mytrainseq.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $genusTaxid;
my $genusName;
my $lastName;
my $nameTest;
my $newLine;
my $j;
my $seq;

#declare array
my @taxon;
my @fasta;
my @line;
my @genusTaxid;
#my @genusName;

#declare hash
my %genusName;

open (TAXON,"<",$ARGV[0]) || die "Error reading mytaxon.txt: $!\n";
@taxon = <TAXON>;
close TAXON;

while($taxon[$i]) {
	$line = $taxon[$i];
	chomp $line;

	if ($line =~ /\*genus/) {
		@line = split(/\*/,$line);
#		$genusTaxid = $line[0];
#		push(@genusTaxid,$genusTaxid);
		$genusName = $line[1];
		$genusName{$genusName} = 1;
	}
	$i++;
	@line=();
	$genusName=();
}
$i=0;

open (FASTA,"<",$ARGV[1]) || die "Error reading mytrainseq.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

open (OUT,">>","mytrainseq_nogenus.fasta") || die "Error writing mytrainseq_nogenus.fasta: $!\n";

while($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/;/,$line);
		$lastName = pop(@line);
		$nameTest = $genusName{$lastName};

		if ($nameTest == 1) {
			$newLine = join(';',@line);
			print OUT $newLine."\n";
			$j=$i+1;
			$seq = $fasta[$j];
			chomp $seq;
			print OUT $seq."\n";
		}
	}
	$i++;
	@line=();
	$lastName=();
	$nameTest=();
	$newLine=();
	$j=();
	$seq=();
}
$i=0;

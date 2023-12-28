#!/usr/bin/perl
#Modified to remove genus, family, order, and class refs so that trained rank eq phylum
#Modified to remove genus, family, and order refs so that trained rank eq class
#Modified to remove genus and family refs so that set is trained to order rank
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
my $familyTaxid;
my $familyName;
my $orderTaxid;
my $orderName;
my $classTaxid;
my $className;

#declare array
my @taxon;
my @fasta;
my @line;

#declare hash
my %genusName;
my %familyName;
my %orderName;
my %className;

open (TAXON,"<",$ARGV[0]) || die "Error reading mytaxon.txt: $!\n";
@taxon = <TAXON>;
close TAXON;

while($taxon[$i]) {
	$line = $taxon[$i];
	chomp $line;

	if ($line =~ /\*genus/) {
		@line = split(/\*/,$line);
		$genusName = $line[1];
		$genusName{$genusName} = 1;
	}
	elsif ($line =~ /\*family/) {
		@line = split(/\*/,$line);
		$familyName = $line[1];
		$familyName{$familyName} = 1;
	}
	elsif ($line =~ /\*order/) {
		@line = split(/\*/,$line);
		$orderName = $line[1];
		$orderName{$orderName} = 1;
	}
	elsif ($line =~ /\*class/) {
		@line = split(/\*/,$line);
		$className = $line[1];
		$className{$className} = 1;
	}

	$i++;
	@line=();
	$genusName=();
	$familyName=();
	$orderName=();
	$className=();
}
$i=0;

open (FASTA,"<",$ARGV[1]) || die "Error reading mytrainseq.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

open (OUT,">>","mytrainseq_nogenus_nofamily_noorder_noclass.fasta") || die "Error writing mytrainseq_nogenus_nofamily_noorder_noclass.fasta: $!\n";

while($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/;/,$line);
		$lastName = pop(@line); #remove genus
		$lastName = pop(@line); #remove family
		$lastName = pop(@line); #remove order
		$lastName = pop(@line); #remove class
		$nameTest = $className{$lastName};

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

#!/usr/bin/perl
#May 13, 2013 edit to target family rank as lowest rank
#March 22, 2013 by Terri Porter
#Script to create a properly formatted fasta file to train the Ribosomal Database Project Naive Bayesian Classifier
#usage perl make_NBC_fasta.plx filtered_taxids.txt.uniq gb_seq.map.uniq taxid.parsed gb_taxid.map.uniq

use strict;
use warnings;

#declare variables
my $i=0;
my $line;
my $tax;
my $seq;
my $lineage;
my $gb;
my $gbLine;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $original;
my $new;
my $j=0;

#declare arrays
my @tax; #filtered list of taxids
my @seq;
my @lineage;
my @line;
my @taxmap; #gb-taxid mapping
my @gbLine;

#declare hashes
my %seq; #indexed by gb
my %lineage; #indexed by taxid
my %taxmap; #indexed by tax

open (TAXIDS, "<", $ARGV[0]) || die "Error cannot open filtered_taxids.txt: $!\n";
@tax = <TAXIDS>;
close TAXIDS;

open (SEQ, "<", $ARGV[1]) || die "Error cannot open gb_seq.map.uniq: $!\n";
@seq = <SEQ>;
close SEQ;

open (LINEAGE, "<", $ARGV[2]) || die "Error cannot open taxid.parsed: $!\n";
@lineage = <LINEAGE>;
close LINEAGE;

open (TAXMAP, "<", $ARGV[3]) || die "Error cannot open gb_taxid.map.uniq: $!\n";
@taxmap = <TAXMAP>;
close TAXMAP;


#hash sequences
while ($seq[$i]) {
	$line = $seq[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/, $line);
		$gb = $line[0];
		$seq = $line[1];
		$seq{$gb} = $seq;
	}
	$i++;
	$line=();
	@line=();
	$gb=();
	$seq=();
}
$i=0;

#hash lineage
while ($lineage[$i]) {
	$line = $lineage[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/, $line);
		$tax = $line[0];
		$kingdom = $line[1];
		$phylum = $line[2];
		$class = $line[3];
		$order = $line[4];
		if ($order eq 'undef') {
			$order = $order."_".$class;
		}
		$family = $line[5];
		if ($family eq 'undef') {
			$family = $family."_".$order;
		}
#		$genus = $line[6];
#		if ($genus eq 'undef') {
#			$i++;
#			$line=();
#			@line=();
#			$tax=();
#			$kingdom=();
#			$phylum=();
#			$class=();
#			$order=();
#			$family=();
#			$genus='OMITTHISONE';
#		}
#		if ($genus eq 'Eutrapela') { ### can be in family Geometridae or Tenebrionidae
#			$genus = $genus."_".$family;
#		}
#		$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus;
		$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family;

		$lineage{$tax} = $lineage;
	}
	$i++;
	$line=();
	@line=();
	$tax=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
#	$genus=();
	$lineage=();
}
$i=0;

#hash gb_taxid.map.uniq ### one taxid can be associated with many different gb ###
while ($taxmap[$i]) {
	$line = $taxmap[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/, $line);
		$gb = $line[0];
		$tax = $line[1];
		if (exists $taxmap{$tax}) {
			$original = $taxmap{$tax};
			$new = $original."|".$gb;
			$taxmap{$tax} = $new;
		}
		else {
			$taxmap{$tax} = $gb;
		}
		$original=();
		$new=();
	}
	$i++;
	$line=();
	@line=();
	$gb=();
	$tax=();
}
$i=0;

open (OUT, ">>", "testNBC.fasta") || die "Error cannot open outfile: $!\n";

while ($tax[$i]) {
	$tax = $tax[$i];
	chomp $tax;

	if (exists $taxmap{$tax}) {
		$gbLine = $taxmap{$tax};
		@gbLine = split(/\|/, $gbLine);

		while ($gbLine[$j]) {
			$gb = $gbLine[$j];

			if (exists $seq{$gb}) {
				$seq = $seq{$gb};
				$seq =~ tr/[A,C,T,G]/[a,c,t,g]/; #to match RDP NBC sample file
	
				if (exists $lineage{$tax}) {
					$lineage = $lineage{$tax};
					if ($lineage =~ /OMITTHISONE/) {
						$j++;
						next;
					}
					else {
						print OUT ">$gb $lineage\n$seq\n";
						delete $seq{$gb};
					}
				}
				else {
					print "Taxid $tax missing from taxid.parsed\n";
				}
			}
			else {
#				print "GB $gb missing from gb_seq.map.uniq\n"; #no need to really print this
			}	
			$j++;
			$gb=();
			$seq=();
			$lineage=();
		}
		$j=0;
		@gbLine=();
	}
	else {
		print "Taxid $tax missing from gb_taxid.map.uniq\n";
	}
	$i++;
	$tax=();
	$gbLine=();
}
$i=0;

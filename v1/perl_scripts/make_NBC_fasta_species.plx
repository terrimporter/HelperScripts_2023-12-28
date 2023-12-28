#!/usr/bin/perl
#Jan. 5, 2017 edit to retain species names for CO1v3_species training set
#Aug. 11, 2016 resolve warnings re: uninitialized rank variables, remaining warnings now refer to taxonomic classifications to poorly identified species so okay to exclude
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
my $species;
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

open (SEQ, "<", $ARGV[1]) || die "Error cannot open gb_seq.uniq: $!\n";
@seq = <SEQ>;
close SEQ;

open (LINEAGE, "<", $ARGV[2]) || die "Error cannot open taxid.parsed: $!\n";
@lineage = <LINEAGE>;
close LINEAGE;

open (TAXMAP, "<", $ARGV[3]) || die "Error cannot open gb_taxid.uniq: $!\n";
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

		if (! defined $tax) {
			print "Taxid not found in taxid.parsed file\n";
		}

		$kingdom = $line[1];

		if (! defined $kingdom) {
			$kingdom = 'undef';
		}
		
		$phylum = $line[2];

		if (! defined $phylum) {
			$phylum = 'undef';
		}
		
		$class = $line[3];

		if (! defined $class) {
			$class = 'undef';
		}

		$order = $line[4];
		
		if (length $order) {

			if ($order eq 'undef') {
				$order = $order."_".$class;
			}
			if ($order eq 'Plecoptera') { ### can be genus of moths or order of stoneflies
				$order = $order."_".$class;
			}
			if ($order eq 'Protura') { ### can be class and order
				$order = $order."_".$class;
			}
			if ($order eq 'Diplura') { ### can be class and order
				$order = $order."_".$class;
			}
		}
		else {
#			print "Problem with order for $tax\n";
			$order = 'undef';
		}

		$family = $line[5];

		if (length $family) {
			if ($family eq 'undef') {
				$family = $family."_".$order;
			}
		}
		else {
#			print "Problem with family for $tax\n";
			$family = 'undef';
		}

		$genus = $line[6];

		if (length $genus) {
			if ($genus eq 'undef') { #this is leading to problems
				print "Genus missing from $tax so rank var changed to undef\n";
				$lineage{$tax} = "undef;undef;undef;undef;undef;OMITTHISONE";
			}
			elsif ($genus eq 'Eutrapela') { ### can be in family Geometridae or Tenebrionidae
				$genus = $genus."_".$family;
				$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus;
				$lineage{$tax} = $lineage;
			}
			elsif ($genus eq 'Plecoptera') { ### can be genus of moths or order of stoneflies
				$genus = $genus."_".$family;
				$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus;
				$lineage{$tax} = $lineage;
			}
			else {
				$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus;
				$lineage{$tax} = $lineage;
			}
		}
		else {
			print "Problem with genus for $tax\n";
			$genus='OMITTHISONE';
			print "Genus missing from $tax so rank var changed to undef (2)\n";
			$lineage{$tax} = "undef;undef;undef;undef;undef;OMITTHISONE";
			$lineage{$tax} = $lineage;
		}

		$species = $line[7];
		$species =~ s/ /_/g; #replace spaces with underscore

		if (length $species) {
			if ($species eq 'undef') {
				print "Species is missing from $tax so rank var changed to undef\n";
				$lineage{$tax} = "undef;undef;undef;undef;undef;undef;OMITTHISONE";
			}
			else {
				$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
				$lineage{$tax} = $lineage;
			}
		}
		else {
			print "Problem with species for $tax\n";
			$species='OMITTHISONE';
			print "Species missing from $tax so rank var changed to undef (3)\n";
			$lineage{$tax} = "undef;undef;undef;undef;undef;undef;OMITTHISONE";
			$lineage{$tax} = $lineage;
		}
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
	$genus=();
	$species=();
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

open (OUT, ">>", "testNBC_species.fasta") || die "Error cannot open outfile: $!\n";

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
					print "Taxid $tax missing from lineage hash\n";
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

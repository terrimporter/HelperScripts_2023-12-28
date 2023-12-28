#!/usr/bin/perl
#Feb. 9/23 edit to handle blank superkingdom
#Apr. 10, 2018 edit to exclude GBs that contain non-nucleotide characters (A,C,G,T); handle 'cellularOrganisms'
#Jan. 5, 2017 edit to retain species rank
#OCt. 25, 2016 edit to add superkingdom rank
#Aug. 11, 2016 resolve warnings re: uninitialized rank variables, remaining warnings now refer to taxonomic classifications to poorly identified species so okay to exclude
#March 22, 2013 by Terri Porter
#Script to create a properly formatted fasta file to train the Ribosomal Database Project Naive Bayesian Classifier
#USAGE perl make_NBC_fasta.plx filtered_taxids.txt.uniq gb_seq.map.uniq taxid.parsed gb_taxid.map.uniq gb_poorqual.txt

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
my $cellularOrganisms;
my $superkingdom;
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
my @gb_poorqual; #filter these poor quality gb's out

#declare hashes
my %seq; #indexed by gb
my %lineage; #indexed by taxid
my %taxmap; #indexed by tax
my %gb_poorqual; #indexed by gb

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

open (POORQUAL, "<", $ARGV[4]) || die "Error cannot open gb_poorqual.txt: $!\n";
@gb_poorqual = <POORQUAL>;
close POORQUAL;

#hash gb_poorqual
while ($gb_poorqual[$i]) {
	$line = $gb_poorqual[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		$line =~ s/ /_/g; #turn any spaces into underscores
		$gb_poorqual{$line} = 1;
	}
	$i++;
	$line=();
}
$i=0;

#hash sequences, but omit the ones that are in gb_poorqual
while ($seq[$i]) {
	$line = $seq[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/, $line);
		$gb = $line[0];
		$gb =~ s/ /_/g; #turn any spaces into underscores
	
		if (exists $gb_poorqual{$gb}) {
			$i++;
			next;
		}
		else {
			$seq = $line[1];
			$seq{$gb} = $seq;
		}
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

		if ($tax eq '') {
			print "Taxid not found in taxid.parsed file\n";
			$i++;
			$line=();
			@line=();
			$tax=();
			$cellularOrganisms=();
			$superkingdom=();
			$kingdom=();
			$phylum=();
			$class=();
			$order=();
			$family=();
			$genus=();
			$species=();
			$lineage=();

			next;

		}
		
		$cellularOrganisms = $line[1];

		$superkingdom = $line[2];

		if ($superkingdom eq '') {
			$lineage{$tax} = "undef;undef;undef;undef;undef;undef;undef;undef;OMITTHISONE";
			print "Superkingdom missing from $tax so flag to OMIT\n";
			$i++;
			$line=();
			@line=();
			$tax=();
			$cellularOrganisms=();
			$superkingdom=();
			$kingdom=();
			$phylum=();
			$class=();
			$order=();
			$family=();
			$genus=();
			$species=();
			$lineage=();

			next;
		}
		
		$kingdom = $line[3];

		if ($kingdom eq '') {
			$kingdom = 'undef_'.$superkingdom;
		}
		
		$phylum = $line[4];

		if ($phylum eq '') {
			$phylum = 'undef_'.$kingdom;
		}
			
		$class = $line[5];

		if ($class eq '') {
			$class = 'undef_'.$phylum;
		}
			
		$order = $line[6];
		
		if ($order eq '') {
			$order = 'undef'.'_'.$class;
		}

		$family = $line[7];

		if ($family eq '') {
			$family = 'undef'.'_'.$order;
		}

		$genus = $line[8];

		if ($genus eq '') { #this is leading to problems
			$lineage{$tax} = "undef;undef;undef;undef;undef;undef;undef;undef;OMITTHISONE";
			print "Genus missing from $tax so flag to OMIT\n";
			$i++;
			$line=();
			@line=();
			$tax=();
			$cellularOrganisms=();
			$superkingdom=();
			$kingdom=();
			$phylum=();
			$class=();
			$order=();
			$family=();
			$genus=();
			$species=();
			$lineage=();

			next;

		}
		
		$species = $line[9];

		if ($species eq '') {
			print "Species missing from $tax so rank var changed to undef\n";
			$lineage{$tax} = "undef;undef;undef;undef;undef;undef;undef;undef;OMITTHISONE";
		}
		else {
			$species =~ s/\s{1}/_/g; # replace space with underscore
			$lineage = $cellularOrganisms.";".$superkingdom.";".$kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
			$lineage{$tax} = $lineage;
		}

	}
		
	$i++;
	$line=();
	@line=();
	$tax=();
	$cellularOrganisms=();
	$superkingdom=();
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
		$gb =~ s/ /_/g; #turn spaces into underscore
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
			$gb =~ s/ /_/g; #replace any spaces with underscores

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
						$lineage =~ s/\s{1}/_/g; #replace any spaces with underscores
						print OUT ">$gb $lineage\n$seq\n";
						delete $seq{$gb};
					}
				}
				else {
					print "Taxid $tax missing from lineage hash\n";
				}
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

#!/usr/bin/perl
#Jan. 5, 2017 edit to retain species rank
#OCt. 25, 2016 edit to add superkingdom rank
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

		$superkingdom = $line[1];

		if (! defined $superkingdom) {
			$superkingdom = 'undef';
		}
		
		$kingdom = $line[2];

		if (length $kingdom) {
			if ($kingdom eq 'undef') {
				$kingdom = $kingdom.'_'.$superkingdom;
			}
		}
		else {
			$kingdom = 'undef';
		}
		
		$phylum = $line[3];

		if (length $phylum) {
			if ($phylum eq 'undef') {
				$phylum = $phylum.'_'.$kingdom;
			}
		}
		else {
			$phylum = 'undef';
		}
		
		$class = $line[4];

		if (length $class) {

			if ($class eq 'undef') {
				$class = $class.'_'.$phylum;
			}
		}
		else {
			$class = 'undef';
		}

		$order = $line[5];
		
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
			$order = 'undef';
		}

		$family = $line[6];

		if (length $family) {
			if ($family eq 'undef') {
				$family = $family."_".$order;
			}
		}
		else {
			$family = 'undef';
		}

		$genus = $line[7];

		if (length $genus) {
			if ($genus eq 'undef') { #this is leading to problems
				$genus = $genus."_".$family;
			}
#				print "Genus missing from $tax so rank var changed to undef\n";
#				$lineage{$tax} = "undef;undef;undef;undef;undef;undef;OMITTHISONE";
#			}
			elsif ($genus eq 'Eutrapela') { ### can be in family Geometridae or Tenebrionidae
				$genus = $genus."_".$family;
			}
			elsif ($genus eq 'Plecoptera') { ### can be genus of moths or order of stoneflies
				$genus = $genus."_".$family;
			}
			elsif ($genus eq 'Alaria') { ### can be a genus of stramenopiles or platyhelminths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Automolus') { ### genus of birds and genus of beetles
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Achlya') { ### genus of oomycetes and lepidoptera
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Ozophora') { ### genus of insecta and red algae
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Acrotylus') { ### genus of red algae and insects
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Lessonia') { ### genus of birds and brown algae
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Chondracanthus') { ### genus of red algae and copepods
				$genus = $genus.'_'.$family;
			}
		}

		$species = $line[8];
		$species =~ s/\s+/_/g; #change spaces to underscore just in case

		if (length $species) {
			if ($species eq 'undef') {
				print "Species missing from $tax so rank var changed to undef\n";
				$lineage{$tax} = "undef;undef;undef;undef;undef;undef;undef;OMITTHISONE";
			}
#			elsif ($species eq 'Bufotes_raddei') { ### compensate for taxonomy issues, Bufotes has no rank
#				$genus = 'Bufotes';
#				$lineage = $superkingdom.";".$kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
#				print "TEST: $lineage\n"; #test
#				$lineage{$tax} = $lineage;
#			}
			else {
				$lineage = $superkingdom.";".$kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
				$lineage{$tax} = $lineage;
			}
		}

		else {
			print "Problem with species for $tax\n";
			$genus='OMITTHISONE';
			print "Species missing from $tax so rank var changed to undef (3)\n";
			$lineage{$tax} = "undef;undef;undef;undef;undef;undef;undef;OMITTHISONE";
			$lineage{$tax} = $lineage;
		}
	}
		
	$i++;
	$line=();
	@line=();
	$tax=();
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

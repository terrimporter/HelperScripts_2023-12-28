#!/usr/bin/perl
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

		if (! defined $tax) {
			print "Taxid not found in taxid.parsed file\n";
		}
		
		$cellularOrganisms = $line[1];

		$superkingdom = $line[2];

		if (! defined $superkingdom) {
			$superkingdom = 'undef';
		}
		
		$kingdom = $line[3];

		if (length $kingdom) {
			if ($kingdom eq 'undef') {
				$kingdom = $kingdom.'_'.$superkingdom;
			}
		}
		else {
			$kingdom = 'undef_'.$superkingdom;
			print "got here\n\n";
		}
		
		$phylum = $line[4];

		if (length $phylum) {
			if ($phylum eq 'undef') {
				$phylum = $phylum.'_'.$kingdom;
			}
#			if ($phylum eq 'Acanthocephala') { # phylum of worms, or genus of insects
#				$phylum = $phylum."_".$kingdom;
#			}
#			if ($phylum eq 'Rhodophyta') { # I assumed everything in BOLD was Metazoa but not
#				$kingdom = 'undef_'.$superkingdom;
#			}
#			if ($phylum eq 'Ctenophora') { ### can be a phylum of Metazoa or genus of arthropods
#				$phylum = $phylum.'_'.$kingdom;
#			}
#			if ($phylum eq 'Actinobacteria') { ### can be phylum and class of Bacteria
#				$phylum = $phylum.'_'.$kingdom;
#			}
#			if ($phylum eq 'Gemmatimonadetes') { ### phylum and class of bacetria
#				$phylum = $phylum.'_'.$kingdom;
#			}
#			if ($phylum eq 'Deferribacteres') { ### phylum and class of bacetria
#				$phylum = $phylum.'_'.$kingdom;
#			}

		}
		else {
			$phylum = 'undef'.'_'.$kingdom;
		}
		
		$class = $line[5];

		if (length $class) {

			if ($class eq 'undef') {
				$class = $class.'_'.$phylum;
			}
#			if ($class eq 'Deferribacteres') { # phylum and class of bacteria
#				$class = $class.'_'.$phylum;
#			}

		}
		else {
			$class = 'undef'.'_'.$phylum;
		}

		$order = $line[6];
		
		if (length $order) {

			if ($order eq 'undef') {
				$order = $order."_".$class;
			}
#			if ($order eq 'Plecoptera') { ### can be genus of moths or order of stoneflies
#				$order = $order."_".$class;
#			}
#			if ($order eq 'Protura') { ### can be class and order
#				$order = $order."_".$class;
#			}
#			if ($order eq 'Diplura') { ### can be class and order
#				$order = $order."_".$class;
#			}
			if ($order eq 'Pristiformes/Rhiniformes group') { ### fix / and space
				$order = 'Pristiformes_Rhiniformes_group';
			}
#			if ($order eq 'Parachela') { ### order of tardigrades and genus of fish
#				$order = $order.'_'.$class;
#			}
			if ($order eq 'Bacteroidetes Order II. Incertae sedis') { ### fix spacing
				$order = 'Bacteroidetes_Order_II_Incertae_sedis';
			}
#			if ($order eq 'Pygophora') { # genus flies, order crustaceans
#				$order = $order.'_'.$class;
#			}

		}
		else {
			$order = 'undef'.'_'.$class;
		}

		$family = $line[7];

		if (length $family) {
			if ($family eq 'undef') {
				$family = $family."_".$order;
			}
#			if ($family eq 'Cepheidae') {
#				$family = $family."_".$order;
#			}
			if ($family eq 'Chilodontidae') { ### family fish and molluscs
				$family = $family.'_'.$order;
			}
		}
		else {
			$family = 'undef'.'_'.$order;
		}

		$genus = $line[8];

		if (length $genus) {
			if ($genus eq 'undef') { #this is leading to problems
				$genus = $genus."_".$family;
			}
			elsif ($genus eq 'Caviria') { # genus in Lymantriidae in NCBI but Erebidae in BOLD
				$family = 'Lymantriidae';
			}
			elsif ($genus eq 'Antachara') { # genus in Noctuidae in NCBI but Erebidae in BOLD
				$family = 'Noctuidae';
			}
			elsif ($genus eq 'Eisenia') { # genus brown algae and segmented worms
				$genus = $genus."_".$family;
			}
			elsif ($genus eq 'Ricinus') { # genus lice and eudicots
				$genus = $genus."_".$family;
			}
			elsif ($genus eq 'Lactarius') { # genus bony fishes and basidiomycetes
				$genus = $genus."_".$family;
			}

			elsif ($genus eq 'Eutrapela') { ### can be in family Geometridae or Tenebrionidae
				$genus = $genus."_".$family;
			}
#			elsif ($genus eq 'Plecoptera') { ### can be genus of moths or order of stoneflies
#				$genus = $genus."_".$family;
#			}
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
			elsif ($genus eq 'Elachista') { ### genus of stramenopiles and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Bostrychia') { ### genus of stoneflies and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Roya') { ### genus of green algae and molluscs
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Rhizophagus') { ### genus of beetles and fungi
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Nemastoma') { ### genus of arthropod and red alga
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Ariopsis') { ### genus of plants and catfish
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Pontogeneia') { ### genus of fungi and amphipods
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Chondria') { ### genus red algae and beetles
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Ptilophora') { ### genus of red algae and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Karenia') { ### genus of alveolata and cicada
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Lobophora') { ### genus of stramenopiles and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Candida') { ### two clades of Candida yeasts
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Dacrydium') { ### genus of conifer and bivalve
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Trigonostomum') { ### genus of flatworms and beetles
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Nectria') { ### genus of fungi and echinodermata
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Dracunculus') { ### genus of plants and nematodes
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Ganonema') { ### genus of red algae and caddisflies
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Bouchetia') { ### genus of plant and sea snail
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Agardhiella') { ### genus of red algae and land snails
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Solieria') { ### genus of red algae and flies
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Diplura') { ### genus of arthropods and stramenopiles
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Zonaria') { ### genus of stramenopiles and gastropods
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Vestia') { ### genus of plants and snails
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Hillia') { ### genus of plants and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Olea') { ### genus of plants and sea slugs
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Grania') { ### genus of annelid worms and red algae
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Atractomorpha') { ### genus of plants and insects
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Bartramia') { ### genus of mosses and birds
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Heterococcus') { ### genus of stramenopiles and arthropods
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Tephrosia') { ### genus of plants and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Rondeletia') { ### genus of plants and fish
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Turbinaria') { ### genus of cnidarians and stramenopiles
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Mastophora') { ### genus of arthropods and red algae
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Drymonia') { ### genus of plants and moths
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Planococcus') { ### genus of bacteria and insects
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Paracoccus') { ### genus of bacteria and mealy bugs
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Rhodococcus') { ### genus of bacteria and scale insect
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Bacillus') { ### genus of bacillus andstick insect
				$genus = $genus.'_'.$family;
			}
			elsif ($genus eq 'Darwinella') { ### genus of beetles and sponges
				$genus = $genus.'_'.$family;
			}
#			elsif ($genus eq 'Pygophora') { # genus flies, order crustaceans
#				$genus = $genus.'_'.$family;
#			}
			elsif ($genus eq 'ApallatesÂ') {
				$genus = 'Apallates'; # handle weird hidden character introduced in BOLD taxon
			}
		}

		$species = $line[9];
#		$species =~ s/\s+/_/g; #change spaces to underscore just in case

		if (length $species) {
			if ($species eq 'undef') {
				print "Species missing from $tax so rank var changed to undef\n";
				$lineage{$tax} = "r__undef;sk__undef;k__undef;p__undef;c__undef;o__undef;f__undef;g__undef;OMITTHISONE";
			}
			elsif ($species =~ /Eutrapela clemataria/) {
				$genus = "Eutrapela_Geometridae";
			}
			elsif ($species =~ /Pelophila borealis/) { # Misannotated as flatworm in GenBank, correct in BOLD, probably problem when depositing to GenBank from BOLD
#				print "Found it\n\n";
				$phylum = 'Arthropoda';
				$class = 'Insecta';
				$order = 'Coleoptera';
				$family = 'Carabidae';
				$genus = 'Pelophila';
				$lineage = "r__$cellularOrganisms;sk__$superkingdom;k__$kingdom;p__$phylum;c__$class;o__$order;f__$family;g__$genus;s__$species";
				$lineage{$tax} = $lineage;
			}
			elsif ($species =~ /coxendix$/) { #handle weird hidden character introduced in BOLD taxon
				$species = 'Apallates coxendix';
				$species =~ s/\s{1}/_/g; # replace any space with underscore
				$lineage = "r__$cellularOrganisms;sk__$superkingdom;k__$kingdom;p__$phylum;c__$class;o__$order;f__$family;g__$genus;s__$species";
				$lineage{$tax} = $lineage;
			}
			else {
				$species =~ s/\s{1}/_/g; # replace space with underscore
				$lineage = "r__$cellularOrganisms;sk__$superkingdom;k__$kingdom;p__$phylum;c__$class;o__$order;f__$family;g__$genus;s__$species";
				$lineage{$tax} = $lineage;
			}
		}

		else {
			print "Problem with species for $tax\n";
			$genus='OMITTHISONE';
			print "Species missing from $tax so rank var changed to undef (3)\n";
			$lineage = "r__undef;sk__undef;k__undef;p__undef;c__undef;o__undef;f__undef;g__undef;OMITTHISONE";
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

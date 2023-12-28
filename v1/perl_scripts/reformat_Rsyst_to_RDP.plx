#!/usr/bin/perl
# Reformat Rsyst:diatom data from parse_Rsyst_diatom.R for the RDP classifier
# USAGE perl reformat_Rsyst_to_RDP.plx SSU.csv

use strict;
use warnings;
use Data::Dumper;

#var
my $i=0;
my $line;
my $id;
my $seq;
my $cellularOrganisms="";
my $domain="";
my $kingdom="";
my $subkingdom="";
my $phylum="";
my $class="";
my $order="";
my $family="";
my $genus="";
my $species="";
my $lineage;
my $acc;
my $fulltaxon;
my $taxon;
my $outfile1 = "testNBC.fasta";
my $outfile2 = "testNBC.taxonomy";
my $prevtaxon;
my $rank;
my $target;
my $termcounter=0;
my $prevtermcounter;
my $lineageindex=0;
my $taxline;
my $lasttaxon;
my $key;
my $value;

#arrays
my @rsyst;
my @fasta;
my @line;
my @lineage;
my @taxon;
my @targets = ("cellularOrganisms","domain","kingdom","subkingdom","phylum","class","order","family","genus","species");
my @split;
my @value;


#hashes
my %taxline; #key =taxon , value = taxonomy line
my %observed; #key = termcounter, value = taxonomyline
my %sorted; #key = termcounter, value = value

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@rsyst = <IN>;
close IN;

# read in the Rsyst diatom csv file
# hash the taxonomy
# reformat to a FASTA file for RDP classifier

open (OUT, ">", $outfile1) || die "Error cannot open fasta outfile1: $!\n";

#hash the taxonomy file 
while ($rsyst[$i]) {
	$line = $rsyst[$i];
	chomp $line;

	# skip header line
	if ($i == 0) {
		$i++;
		next;
	}
 
	@line = split(/,/,$line);
	$id = $line[0];
	$seq = $line[1];
	$domain = $line[2];
	$kingdom = $line[3];
	$subkingdom = $line[4];
	$phylum = $line[5];
	$class = $line[6];
	$order = $line[7];
	$family = $line[8];
	$genus = $line[9];
	$species = $line[10];

	unless (length $species) {
		print "problem processing id $id\n";
	}

	$species =~ s/ /_/g;
	$species =~ s/\.//g;
	$species =~ s/\(//g;
	$species =~ s/\)//g;

	$lineage = "cellularOrganisms;".$domain.";".$kingdom.";".$subkingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;

	print OUT ">$id $lineage\n$seq\n";

	$i++;
	$domain="";
	$kingdom="";
	$subkingdom="";
	$phylum="";
	$class="";
	$order="";
	$family="";
	$genus="";
	$species="";
}
$i=0;
close OUT;

open (IN2, "<", $outfile1) || die "Cannot open testNBC.fasta: $!\n";
@fasta = <IN2>;
close IN2;

# parse through fasta file to hash taxonomy
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		parse_header();

		foreach my $target (@targets) {
			if (exists $observed{$target}) {
				$taxon = $observed{$target};

				unless (exists $taxline{$taxon}) {
					$termcounter++;
					if ($target eq "cellularOrganisms") {
						$prevtermcounter=$termcounter-1;
						$taxline{$taxon} = $termcounter."*".$taxon."*".$prevtermcounter."*".$lineageindex."*".$target;
					}
					else {
						$taxline = $taxline{$lasttaxon};
						@split = split(/\*/,$taxline);
						$prevtermcounter = $split[0];
						$taxline{$taxon} = $termcounter."*".$taxon."*".$prevtermcounter."*".$lineageindex."*".$target;
					}
				}
				$lasttaxon = $taxon;
				$lineageindex++;
			}
		}
		$lineageindex=0;
	}
	$i++;
	%observed=();
}
$i=0;


#sort the hash by termcounter into a new hash indexed by termcounter
while ( ($key, $value) = each (%taxline) ) {
	@value = split(/\*/, $value);
	$termcounter = $value[0];
	$sorted{$termcounter} = $value;
}
$key=();
$value=();

open (OUT2, ">>", $outfile2) || die "Cannot open taxonomy outfile2: $!\n";

foreach $key (sort {$a <=> $b} keys %sorted) {
	$value = $sorted{$key};
	print OUT2 $value."\n";
}
close OUT2;

##########

sub parse_header {
	
	$line =~ s/^>//g;
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);

	@lineage = split(/;/,$lineage);
	$cellularOrganisms = $lineage[0];
	$domain = $lineage[1];
	$kingdom = $lineage[2];
	$subkingdom = $lineage[3];
	$phylum = $lineage[4];
	$class = $lineage[5];
	$order = $lineage[6];
	$family = $lineage[7];
	$genus = $lineage[8];
	$species = $lineage[9];

	$observed{"cellularOrganisms"} = $cellularOrganisms;
	$observed{"domain"} = $domain;
	$observed{"kingdom"} = $kingdom;
	$observed{"subkingdom"} = $subkingdom;
	$observed{"phylum"} = $phylum;
	$observed{"class"} = $class;
	$observed{"order"} = $order;
	$observed{"family"} = $family;
	$observed{"genus"} = $genus;
	$observed{"species"} = $species;

}

##########



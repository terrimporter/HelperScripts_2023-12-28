#!/usr/bin/perl
# Teresita M. Porter, Oct. 9, 2021
# Edit unite_derep.fasta.strict from UNITE to handle unidentified
# USAGE perl make_NBC_fasta_root_species.plx unite_derep.fasta.strict

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $start_root;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $outfile = "testNBC.fasta";
my $newline;
my $lineage;
my $remove;
my $ref;

# declare arrays
my @in;
my @line;

# declare hash
my %species; #key = species, value = lineage
my %genus;
my %family;
my %order;
my %class;
my %phylum;

open (IN, "<", $ARGV[0]) || die "Error cannnot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

# check for duplicate species
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) { # process header
		@line = split(/;/,$line);
		$start_root = $line[0];
		$kingdom = $line[1];
		$phylum = $line[2];
		$class = $line[3];
		$order = $line[4];
		$family = $line[5];
		$genus = $line[6];
		$species = $line[7];


		$lineage = $kingdom;
		if (exists $phylum{$phylum}) {
			$ref = $phylum{$phylum};

			if ($ref ne $lineage) {
				$phylum = $kingdom."_".$phylum;
			}
		}
		else {
			$phylum{$phylum}=$lineage;
		}


		$lineage = $phylum.";".$class;
		if (exists $class{$class}) {
			$ref = $class{$class};

			if ($ref ne $lineage) {
				$class = $phylum."_".$class;
			}
		}
		else {
			$class{$class}=$lineage;
		}


		$lineage = $kingdom.";".$phylum.";".$class;
		if (exists $order{$order}) {
			$ref = $order{$order};

			if ($ref ne $lineage) {
				$order = $class."_".$order;
			}
		}
		else {
			$order{$order}=$lineage;
		}


		$lineage = $kingdom.";".$phylum.";".$class.";".$order;
		if (exists $family{$family}) {
			$ref = $family{$family};

			if ($ref ne $lineage) {
				$family = $order."_".$family;
			}
		}
		else {
			$family{$family}=$lineage;
		}


		$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family;
		if (exists $genus{$genus}) {
			$ref = $genus{$genus};

			if ($ref ne $lineage) {
				$genus = $family."_".$genus;
			}
		}
		else {
			$genus{$genus}=$lineage;
		}
		if ($genus =~ "g__candida") { # conflict with g__Candida
			$genus = $family."_".$genus;
		}


		$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus;
		if (exists $species{$species}) {
			$ref = $species{$species};

			if ($ref ne $lineage) {
				$species = $genus."_".$species;
			}
		}
		else {
			$species{$species}=$lineage;
		}
		if ($species =~ /s__candida/) { # conflict with s__Candida
			$species = $genus."_".$species;
		}

		print OUT $start_root.";".$kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species."\n";

	}
	else { # process seq
		print OUT $line."\n";
	}
	$i++;
}
$i=0;

close OUT;


#!/usr/bin/perl
#Oct.5, 2010 by Terri Porter
#Script to grab accession and Genus and species from original_genbank.fasta to be used to help parse SAP results
#usage $perl grab_acc_genus_species.plx < original_genbank.fasta > original_genbank.map

use strict;
use warnings;

#declare var
my $line;
my $accession;
my $description;
my $species;

#declare array
my @line;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ />/) {
		@line = split(/\|/,$line);
		$accession = $line[3];
		$description = $line[4];
		if ($description =~ /bayanus/) {
			print "$accession\tSaccharomyces\tbayanus\n";
		}
		elsif ($description =~ /cariocanus/) {
			print "$accession\tSaccharomyces\tcariocanus\n";
		}
		elsif ($description =~ /cerevisiae/) {
			print "$accession\tSaccharomyces\tcerevisiae\n";
		}
		elsif ($description =~ /kudriavzevii/) {
			print "$accession\tSaccharomyces\tkudriavzevii\n";
		}
		elsif ($description =~ /mikatae/) {
			print "$accession\tSaccharomyces\tkudriavzevii\n";
		}
		elsif ($description =~ /paradoxus/) {
			print "$accession\tSaccharomyces\tparadoxus\n";
		}
		elsif ($description =~ /pastorianus/) {
			print "$accession\tSaccharomyces\tpastorianus\n";
		}
		elsif ($description =~ /bulderi/) {
			print "$accession\tKazachstania\tbulderi\n";
		}
		elsif ($description =~ /turicensis/) {
			print "$accession\tKazachstania\tturicensis\n";
		}
		elsif ($description =~ /kunashirensis/) {
			print "$accession\tKazachstania\tkunashirensis\n";
		}
		elsif ($description =~ /rosinii/) {
			print "$accession\tKazachstania\trosinii\n";
		}
		elsif ($description =~ /martiniae/) {
			print "$accession\tKazachstania\tmartiniae\n";
		}
		elsif ($description =~ /transvaalensis/) {
			print "$accession\tKazachstania\ttransvaalensis\n";
		}
		elsif ($description =~ /servazzii/) {
			print "$accession\tKazachstania\tservazzii\n";
		}
		elsif ($description =~ /boulardii/) {
			print "$accession\tKazachstania\tboulardii\n";
		}
		elsif ($description =~ /uvarum/) {
			print "$accession\tuncertain classification\n";
		}
		elsif ($description =~ /ellipsoideus/) {
			print "$accession\tSaccharomyces\tcerevisiae (ellipsoideus)\n";
		}
		elsif ($description =~ /chevalieri/) {
			print "$accession\tuncertain classification\n";
		}
		elsif ($description =~ /spencerorum/) {
			print "$accession\tKazachstania\tspencerorum\n";
		}
		elsif ($description =~ /Saccharomyces sp./) {
			print "$accession\tSaccharomyces\tsp.\n";
		}
		elsif ($description =~ /Saccharomyces \S+/){
			$description =~ /Saccharomyces (\S+)/;
			$species = $1;
			print "$accession\tSaccharomyces\t$species\tmisclassified in Genbank\n";
		}
		elsif ($description =~ /S.\w+/) {
			$description =~ /S.(\w+)/;
			$species = $1;
			print "accession\tSaccharomyces\t$species\tmisclassified in Genbank\n";
		}
		else {
			print "$accession\tmissing\tmissing\n";
		}
	}
}


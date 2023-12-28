#!/usr/bin/perl
#July 27, 2018 include problem_species.txt file (from check_for_SIVLA_inconsistencies.plx)
#July 25, 2018 by Teresita M. Porter
#Script to fix FASTA file so only whole taxon ranks are represented (domain, kingdom, phylum, class, order, family, genus, species)
#convert U's to T's, in lower case
#ensure output is in strict FASTA format (no line breaks in seq)
#USAGE perl fix_SILVA_fasta_headers.plx tax_slv_ssu_132.txt SILVA_132_SSUParc_tax_silva.fasta.gz problem_species.txt

use strict;
use warnings;

#var
my $i=0;
my $line;
my $path;
my $rank;
my $acc;
my $lineage;
my $fulltaxon;
my $taxon;
my $prevtaxon="";
my $newlineage;
my $outfile = "Eukaryota.fasta";
my $flag=0; #flag start of seq
my $oldline="";
my $euk=0; #flag to grab just Eukaryota
my $domain;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $prob=0; #flag record to skip

#arrays
my @tax;
my @fasta;
my @line;
my @lineage;
my @taxon;
my @problem;

#hashes
my %map; #key=path, value=rank
my %duplicates; #key=taxon, value=1
my %target; #key=rank, value=1
my %observed; #key = observed rank in SILVA_132_SSUParc_tax_silva.fasta.gz, value = observed taxon name
my %problem; #key = species, value = 1

#assign ranks to %target
%target = ( #target ranks from tax_slv_ssu_132.txt 
		"domain" => 1,
		"kingdom" => 1,
		"phylum" => 1,
		"class" => 1,
		"order" => 1,
		"family" => 1,
		"genus" => 1 );

#handle duplicate taxon names
%duplicates= (
		"marine_group" => 1,
		"undef" => 1,
		"unidentified" => 1,
		"uncultured" => 1,
		"Incertae_Sedis" => 1,
		"Candida" => 1,
		"Rhizophagus" => 1,
		"Navicula" => 1,
		"Karenia" => 1,
		"Nectria" => 1,
		"Lobophora" => 1,
		"Zonaria" => 1,
		"Ctenophora" => 1,
		"Campanella" => 1,
		"Uronema" => 1,
		"Lessonia" => 1,
		"Atractomorpha" => 1,
		"Bostrychia" => 1,
		"Pontogeneia" => 1,
		"Mastophora" => 1,
		"Protura" => 1,
		"Parachela" => 1,
		"Acrotylus" => 1,
		"Plecoptera" => 1,
		"Achlya" => 1,
		"Nemastoma" => 1,
		"Choanoflagellida" => 1,
		"Ptilophora" => 1,
		"Olea" => 1,
		"Pterospora" => 1,
		"Dacrydium" => 1,
		"Beauveria" => 1,
		"Phoma" => 1,
		"Drymonia" => 1,
		"Roya" => 1,
		"Heterococcus" => 1,
		"Actinobacteria" => 1,
		"Rhodococcus" => 1,
		"Bacillus" => 1,
		"Family_XI" => 1,
		"Unknown_Family" => 1,
		"Thermotogae" => 1,
		"Aquificae" => 1,
		"Family_XII" => 1,
		"Deferribacteres" => 1,
		"Chlamydiae" => 1,
		"Nitrospira" => 1,
		"Latescibacteria" => 1,
		"Chrysiogenetes" => 1,
		"Gemmatimonadetes" => 1,
		"Planococcus" => 1,
		"Elusimicrobia" => 1,
		"Paracoccus" => 1);

#read in the SILVA taxonomy file
open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@tax = <IN>;
close IN;

#hash the taxonomy file 
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$path = $line[0];
	$rank = $line[2];

	$path =~ s/\s+/_/g;
	$path =~ s/\(//g;
	$path =~ s/\)//g;
	$path =~ s/\</_/g;
	$path =~ s/\>/_/g;

	$map{$path} = $rank;

	$i++;
}
$i=0;

#read in problem_species file that contains list of species with more than one taxonomic path
open (IN2, "<", $ARGV[2]) || die "Cannot open problem species infile: $!\n";
@problem = <IN2>;
close IN2;

#hash the problem species
while ($problem[$i]) {
	$line = $problem[$i];
	chomp $line;

	$line =~ s/\s+/_/g;
	$line =~ s/\(//g;
	$line =~ s/\)//g;
	$line =~ s/\.//g;
	$line =~ s/\</_/g;
	$line =~ s/\>/_/g;

	$problem{$line} = 1;
	$i++;
}
$i=0;

#read in the SILVA fasta.gz file
@fasta = `zcat $ARGV[1]`;

#create an outfile
open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	#parse FASTA header
	if ($line =~ /^>/) {
		if ($line =~ /Eukaryota;/) { #just process Eukaryota
			$euk=1;
		}
		else {
			$euk=0;
			$i++;
			next;
		}

		if ($flag==1 && $euk == 1 ) { #print the whole seq before parsing next header
			if ($prob == 0) { #only print out previous seq if it isn't a problem species

				print OUT "$oldline\n";
				$oldline="";

				parse_header();
				reformat_lineage();

				if ($prob == 0) {
					print OUT ">$acc $newlineage\n";
				}
			}
			else {
				$prob=0;
				parse_header();
				reformat_lineage();

				if ($prob==0) {
					print OUT ">$acc $newlineage\n";
				}
			}
		}
		elsif ($flag==0 && $euk == 1) {#handles just the first record
			$flag=1; #flag start of seq
			parse_header();
			reformat_lineage();
			if ($prob==0) {
				print OUT ">$acc $newlineage\n";
			}
			else {
				$i++;
				next;
			}
		}
		else {
			$i++;
			next;
		}
	}
	#concatenate sequence across lines
	elsif ($euk == 1 && $prob == 0) {
		$line =~ s/U/T/g; #convert U's to T's
		$line = lc $line; #convert to lower case
		$oldline = $oldline.$line; #append new line to old line
	}
	elsif ($euk == 0) {
		$i++;
		next;
	}
	elsif ($prob==1) {
		$i++;
		next;
	}
	$i++;
	%observed=();
}

#don't forget to print last seq
print OUT "$oldline\n";
$i=0;
close OUT;


##########

sub parse_header {

	$line =~ s/^>//g;
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;

	@lineage = split(/;/,$lineage);	
		foreach $taxon (@lineage) {
			$fulltaxon = $prevtaxon.$taxon.";";
			$taxon = $fulltaxon;
			$taxon =~ s/;$//g;
			@taxon = split(/;/,$taxon);
			$taxon = pop (@taxon); #grab the last taxon off the lineage path

			if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
				$rank = $map{$fulltaxon};
				$observed{$rank} = $taxon;
			}
			else {
				$observed{"species"} = $taxon; 
			}
			$prevtaxon = $fulltaxon;
		}
	$prevtaxon="";
}

##########

sub reformat_lineage {
	if (exists $observed{"domain"}) {
		$domain = $observed{"domain"};
		$newlineage = $domain.";";
	}
	else {
		print "Couldn't find domain";
	}

	if (exists $observed{"kingdom"}) {
		$kingdom = $observed{"kingdom"};
		$newlineage = $newlineage.$kingdom.";";
	}
	else {
		$kingdom = $domain."_undef;";
		$newlineage = $newlineage.$kingdom;
	}

	if (exists $observed{"phylum"}) {
		$phylum = $observed{"phylum"};

		if (exists $duplicates{$phylum}) {# handle duplicate phyla
			$kingdom =~ s/;$//;
			$phylum = $kingdom."_".$phylum;
		}
		$newlineage = $newlineage.$phylum.";";
	}
	else {
		$kingdom =~ s/;$//;
		$phylum = $kingdom."_undef;";
		$newlineage = $newlineage.$phylum;
	}

	if (exists $observed{"class"}) {
		$class = $observed{"class"};

		if (exists $duplicates{$class}) { #handle duplicate classes
			$phylum =~ s/;$//;
			$class = $phylum."_".$class;
		}
		$newlineage = $newlineage.$class.";";
	}
	else {
		$phylum =~ s/;$//;
		$class = $phylum."_undef;";
		$newlineage = $newlineage.$class;
	}

	if (exists $observed{"order"}) {
		$order = $observed{"order"};

		if (exists $duplicates{$order}) { #handle duplicate orders
			$class =~ s/;$//;
			$order = $class."_".$order;
		}
		$newlineage = $newlineage.$order.";";

	}
	else {
		$class =~ s/;$//;
		$order = $class."_undef;";
		$newlineage = $newlineage.$order;
	}

	if (exists $observed{"family"}) {
		$family = $observed{"family"};

		if (exists $duplicates{$family}) {#handle duplicate families
			$order =~ s/;$//;
			$family = $order."_".$family;
		}
		$newlineage = $newlineage.$family.";";
	}
	else {
		$order =~ s/;$//;
		$family = $order."_undef;";
		$newlineage = $newlineage.$family;
	}

	if (exists $observed{"genus"}) {
		$genus = $observed{"genus"};

		if (exists $duplicates{$genus}) { #handle duplicate genera
			$family =~ s/;$//;
			$genus = $family."_".$genus;
		}
		$newlineage = $newlineage.$genus.";";
	}
	else {
		$family =~ s/;$//;
		$genus = $family."_undef;";
		$newlineage = $newlineage.$genus;
	}

	if (exists $observed{"species"}) {
		$species = $observed{"species"};

		if (exists $problem{$species}) { #don't process species with more than one taxonomic path
			$prob=1;
		}

		if (exists $duplicates{$species}) { #handle duplicate species
			$genus =~ s/;$//;
			$species = $genus."_".$species;
		}
		$newlineage = $newlineage.$species;
	}
	else {
		$genus =~ s/;$//;
		$species = $genus."_undef;";
		$newlineage = $newlineage.$species;
	}

	$domain=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();

}

##########

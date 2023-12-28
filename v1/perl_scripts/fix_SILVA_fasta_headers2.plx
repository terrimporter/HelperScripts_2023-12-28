#!/usr/bin/perl
#July 25, 2018 by Teresita M. Porter
#Adapted to target Bacteri & Archaea from SSU Ref 99 file, subsampled by clustering at 90% sequence similarity
#Script to fix FASTA file so only whole taxon ranks are represented (domain, kingdom, phylum, class, order, family, genus, species)
#convert U's to T's, in lower case
#ensure output is in strict FASTA format (no line breaks in seq)
#USAGE perl fix_SILVA_fasta_headers.plx tax_slv_ssu_132.txt silva.centroids.gz

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
my $outfile = "BacteriaArchaea.fasta";
my $flag=0; #flag start of seq
my $oldline="";
my $pro=0; #flag to grab just prokaryotes
my $domain;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;

#arrays
my @tax;
my @fasta;
my @line;
my @lineage;
my @taxon;

#hashes
my %map; #key=path, value=rank
my %target = ( #target ranks from tax_slv_ssu_132.txt 
		"domain" => 1,
		"kingdom" => 1,
		"phylum" => 1,
		"class" => 1,
		"order" => 1,
		"family" => 1,
		"genus" => 1 );
my %observed; #key = observed rank in SILVA_132_SSUParc_tax_silva.fasta.gz, value = observed taxon name

#handle duplicate taxon names
my %duplicates= (
		"marine_group" => 1);

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
#	print $path."\n"; #test

	$map{$path} = $rank;

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
		if ($line =~ /Bacteria;/ || $line =~ /Archaea;/) { #just process prokaryotes
			$pro=1;
		}
		else {
			$pro=0;
			$i++;
			next;
		}

		if ($flag==1 && $pro == 1) { #print the whole seq before parsing next header
			print OUT "$oldline\n";
			$oldline="";

			#parse header
			parse_header();

			#create properly formatted lineage
			reformat_lineage();

			print OUT ">$acc $newlineage\n";
		}
		elsif ($flag==0 && $pro == 1) {
			$flag=1; #flag start of seq
			parse_header();
			reformat_lineage();
			print OUT ">$acc $newlineage\n";
		}
	}
	#concatenate sequence across lines
	elsif ($pro == 1) {
		$line =~ s/U/T/g; #convert U's to T's
		$line = lc $line; #convert to lower case
		$oldline = $oldline.$line; #append new line to old line
	}
	elsif ($pro == 0 ) {
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
	$line =~ s/_/ /; #replace first underscore with space
	print $line."\n";
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods

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
		$newlineage = $newlineage."undef;";
	}

	if (exists $observed{"phylum"}) {
		$phylum = $observed{"phylum"};
		$newlineage = $newlineage.$phylum.";";
	}
	else {
		$newlineage = $newlineage."undef;";
	}

	if (exists $observed{"class"}) {
		$class = $observed{"class"};
		$newlineage = $newlineage.$class.";";
	}
	else {
		$newlineage = $newlineage."undef;";
	}

	if (exists $observed{"order"}) {
		$order = $observed{"order"};
		$newlineage = $newlineage.$order.";";
	}
	else {
		$newlineage = $newlineage."undef;";
	}

	if (exists $observed{"family"}) {
		$family = $observed{"family"};
		$newlineage = $newlineage.$family.";";
	}
	else {
		$newlineage = $newlineage."undef;";
	}

	if (exists $observed{"genus"}) {
		$genus = $observed{"genus"};

		if (exists $duplicates{$genus}) { #handle duplicate genera
			$genus = $family."_".$genus;
		}

		$newlineage = $newlineage.$genus.";";
	}
	else {
		$newlineage = $newlineage."undef;";
	}

	if (exists $observed{"species"}) {
		$species = $observed{"species"};
		$newlineage = $newlineage.$species;
	}
	else {
		$newlineage = $newlineage."undef;";
	}

}

##########

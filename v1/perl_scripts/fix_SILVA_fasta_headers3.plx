#!/usr/bin/perl
# Jan. 15/20 edit to accomodate a mapping file of problematic lineages
#July 27, 2018 include problem_species.txt file (from check_for_SIVLA_inconsistencies.plx)
#July 25, 2018 by Teresita M. Porter
#Script to fix FASTA file so only whole taxon ranks are represented (domain, kingdom, phylum, class, order, family, genus, species)
#convert U's to T's, in lower case
#ensure output is in strict FASTA format (no line breaks in seq)
#USAGE perl fix_SILVA_fasta_headers.plx tax_slv_ssu_138.txt SILVA_138_SSURef_NR99_tax_silva.fasta.gz problem_map.txt

use strict;
use warnings;
use Data::Dumper;

#var
my $i=0;
my $line;
my $path;
my $lineage;
my $abundance;
my $rank;
my $acc;
my $fulltaxon;
my $taxon;
my $prevtaxon="";
my $newlineage;
my $outfile = "Eukaryota.fasta";
my $seq="";
my $euk=0; #flag to grab just Eukaryota
my $domain;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $prob=0; #flag record to skip, 0 = ok, 1 = in %problem, 2 = using mapping from in %fixproblem
my $good_seq_counter=0;
my $salvaged_seq_counter=0;
my $bad_seq_counter=0;
my $bact_seq_counter=0;

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
my %observed; #key = observed rank in SILVA_132_SSUParc_tax_silva.fasta.gz, value = observed taxon name
my %problem; #outerkey = species, innerkey = lineage, value = abundance
my %fixproblem; #key = species, value = lienage

# handle duplicate taxon names
# add to this list when errors pop up during RDP classifier training
%duplicates= (
		"marine_group;" => 1,
		"undef;" => 1,
		"unidentified;" => 1,
		"uncultured;" => 1,
		"Incertae_Sedis;" => 1,
		"Candida;" => 1,
		"Rhizophagus;" => 1,
		"Navicula;" => 1,
		"Karenia;" => 1,
		"Nectria;" => 1,
		"Lobophora;" => 1,
		"Zonaria;" => 1,
		"Ctenophora;" => 1,
		"Campanella;" => 1,
		"Uronema;" => 1,
		"Lessonia;" => 1,
		"Atractomorpha;" => 1,
		"Bostrychia;" => 1,
		"Pontogeneia;" => 1,
		"Mastophora;" => 1,
		"Protura;" => 1,
		"Parachela;" => 1,
		"Acrotylus;" => 1,
		"Plecoptera;" => 1,
		"Achlya;" => 1,
		"Nemastoma;" => 1,
		"Choanoflagellida;" => 1,
		"Ptilophora;" => 1,
		"Olea;" => 1,
		"Pterospora;" => 1,
		"Dacrydium;" => 1,
		"Beauveria;" => 1,
		"Phoma;" => 1,
		"Drymonia;" => 1,
		"Roya;" => 1,
		"Heterococcus;" => 1,
		"Actinobacteria;" => 1,
		"Rhodococcus;" => 1,
		"Bacillus;" => 1,
		"Family_XI;" => 1,
		"Unknown_Family;" => 1,
		"Thermotogae;" => 1,
		"Aquificae;" => 1,
		"Family_XII;" => 1,
		"Deferribacteres;" => 1,
		"Chlamydiae;" => 1,
		"Nitrospira;" => 1,
		"Latescibacteria;" => 1,
		"Chrysiogenetes;" => 1,
		"Gemmatimonadetes;" => 1,
		"Planococcus;" => 1,
		"Elusimicrobia;" => 1,
		"Paracoccus;" => 1);

#read in the SILVA taxonomy file (domain -> genus)
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
	$path =~ s/\.//g;
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

	@line = split(/\t/, $line);
	$species = $line[0];
	$lineage = $line[1];
	$abundance = $line[2];

	$species =~ s/\s+/_/g;
	$species =~ s/\(//g;
	$species =~ s/\)//g;
	$species =~ s/\.//g;
	$species =~ s/\</_/g;
	$species =~ s/\>/_/g;

	$problem{$species}{$lineage} = $abundance;
	$i++;
}
$i=0;

# remove species from %problem if the dominant lineage can be salvaged into new %fixproblem

# loop through hash of hashes
# find inner key with largest value, and indicate if multiple keys tied for largest value
foreach $species (keys %problem) {

	($newlineage, $prob) = key_with_highest_val($problem{$species});

	if ($prob == 0) { # largest value found, else prob == 1 indicates a tie
		$fixproblem{$species} = $newlineage;
		delete($problem{$species});
	}
}

# read in the SILVA fasta.gz file
@fasta = `zcat $ARGV[1]`;

# create an outfile
open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

# parse through SILVA fasta file
# check for higher level duplicate taxa that needs to be reformatted
# check for problem and fixproblem taxa
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	#parse FASTA header
	if ($line =~ /^>/) {
		if (length $seq) { # check if a sequence exists
			if (($euk == 1 && $prob == 0) || ($euk == 1 && $prob == 2)) {
				print OUT "$seq\n";
				$seq="";
			}
		}

		if ($line =~ /Eukaryota;/) { #just process Eukaryota for now
			$euk=1;

			parse_header();
			reformat_lineage();

			if ($prob == 0) {
				print OUT ">$acc $newlineage\n";
				$good_seq_counter++;
			}
			elsif ($prob == 2) {
				parse_header2();
				reformat_lineage2();
				print OUT ">$acc $newlineage\n";
				$salvaged_seq_counter++;
			}
			else { # prob == 1 indicates species has multiple lineages in SILVA taxonomy file, can't choose a best one (tie among choices)
				$bad_seq_counter++;
				$i++;
				next;
			}
		}
		else {
			$euk=0;
			$bact_seq_counter++;
			$i++;
			next;
		}
	}
	#concatenate sequence across lines
	elsif (($euk == 1 && $prob == 0) || ($euk == 1 && $prob == 2)) {
		$line =~ s/U/T/g; #convert U's to T's
		$line = lc $line; #convert to lower case
		$seq = $seq.$line; #append new line to old line
	}
	$i++;
	%observed=();
	$newlineage=();
}
$i=0;

#don't forget to print last seq if appropriate
if (($euk == 1 && $prob == 0) || ($euk == 1 && $prob == 2)) {
	print OUT "$seq\n";
}
close OUT;

# print status of seqs (good, salvaged, problems)
print "good: $good_seq_counter\n";
print "salvaged: $salvaged_seq_counter\n";
print "bad: $bad_seq_counter\n";
print "bacteria: $bact_seq_counter\n";

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

	$species = pop(@lineage);
	$observed{"species"} = $species;
	
	my @species = split('_',$species); # do this to accomodate for genera missing from SILVA tax file ex. Zea
	$genus = shift(@species);
	$observed{"genus"} = $genus;

	print "testing1\t$genus\t$species\n"; #testing

	my $length = scalar(@lineage); # do this to accomodate for >Accession kingdom;genus_species
	if ($length > 1) {
		my $oldgenus = pop(@lineage); 
	}

	foreach $taxon (@lineage) {
		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
		@taxon = split(/;/,$fulltaxon);
		$taxon = pop(@taxon); # use this to add to %observed

		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
			$rank = $map{$fulltaxon};
			$observed{$rank} = $taxon;
		}
		$prevtaxon = $fulltaxon;
	}
	$prevtaxon="";
}

##########

sub reformat_lineage {
	if (exists $observed{"domain"}) {
		$domain = $observed{"domain"};
		$domain = $domain.";";
	}
	else {
		print "Couldn't find domain for accession $acc\n";
		$domain="unknownDomain;";
	}

	if (exists $observed{"kingdom"}) {
		$kingdom = $observed{"kingdom"};
		$kingdom = $kingdom.";";
	}
	else {
		$domain =~ s/;$//g;
		$kingdom = $domain."_undef;";
		$domain = $domain.";";
	}

	if (exists $observed{"phylum"}) {
		$phylum = $observed{"phylum"};
		$phylum = $phylum.";";

		if (exists $duplicates{$phylum}) {# handle duplicate phyla
			$kingdom =~ s/;$//;
			$phylum = $kingdom."_".$phylum;
			$kingdom = $kingdom.";";
		}
	}
	else {
		$kingdom =~ s/;$//;
		$phylum = $kingdom."_undef;";
		$kingdom = $kingdom.";";
	}

	if (exists $observed{"class"}) {
		$class = $observed{"class"};
		$class = $class.";";

		if (exists $duplicates{$class}) { #handle duplicate classes
			$phylum =~ s/;$//;
			$class = $phylum."_".$class;
			$phylum = $phylum.";";
		}
	}
	else {
		$phylum =~ s/;$//;
		$class = $phylum."_undef;";
		$phylum = $phylum.";";
	}

	if (exists $observed{"order"}) {
		$order = $observed{"order"};
		$order = $order.";";

		if (exists $duplicates{$order}) { #handle duplicate orders
			$class =~ s/;$//;
			$order = $class."_".$order;
			$class = $class.";";
		}
	}
	else {
		$class =~ s/;$//;
		$order = $class."_undef;";
		$class = $class.";";
	}

	if (exists $observed{"family"}) {
		$family = $observed{"family"};
		$family = $family.";";

		if (exists $duplicates{$family}) {#handle duplicate families
			$order =~ s/;$//;
			$family = $order."_".$family;
			$order = $order.";";
		}
	}
	else {
		$order =~ s/;$//;
		$family = $order."_undef;";
		$order = $order.";";
	}

	if (exists $observed{"genus"}) {
		$genus = $observed{"genus"};
		$genus = $genus.";";
	}

	if (exists $observed{"species"}) {
		$species = $observed{"species"};

		if (exists $problem{$species}) { #don't process species with more than one taxonomic path
			$prob=1;
		}

		elsif (exists $fixproblem{$species}) {
			$prob=2;
		}
		else {
			$prob=0;
		}

		$newlineage = $domain.$kingdom.$phylum.$class.$order.$family.$genus.$species;
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

sub parse_header2 {

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

	# parse the original lineage
	@lineage = split(/;/,$lineage);

	$species = pop(@lineage);
	$observed{"species"} = $species;

#	my @species = split('_',$species); # do this to accomodate for genera missing from SILVA tax file ex .Zea
#	$genus = shift(@species);
#	$observed{"genus"} = $genus;

#	my $length = scalar(@lineage); # do this to accomodate for >Accession kingdom;genus_species
#	if ($length > 1) {
#		my $oldgenus = pop(@lineage);
#	}

#	foreach $taxon (@lineage) {
#		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
#		@taxon = split(/;/,$fulltaxon);
#		$taxon = pop(@taxon); # use this to add to %observed

#		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
#			$rank = $map{$fulltaxon};
#			$observed{$rank} = $taxon;
#		}
#		$prevtaxon = $fulltaxon;
#	}

#	$prevtaxon="";	
	
	# grab the most abundant lineage from problem species	
#	$species = $observed{"species"};
	$lineage = $fixproblem{$species};

	@lineage = split(/;/, $lineage); # domain -> genus
	
#	$species = pop(@lineage);
#	$observed{"species"} = $species;

	my @species = split('_',$species); # do this to accomodate for genera missing from SILVA tax file ex. Zea
	$genus = shift(@species);
	$observed{"$genus"} = $genus;

	print "testing2\t$genus\t$species\n"; #testing

	my $length = scalar(@lineage); # do this to accomodate for >Accession kingdom;genus_species
	if ($length > 1) {
		my $oldgenus = pop(@lineage);
	}

	# parse the most abundant lineage
	foreach $taxon (@lineage) {
		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
		@taxon = split(/;/,$fulltaxon);
		$taxon = pop(@taxon); # use this to add to %observed

		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
			$rank = $map{$fulltaxon};
			$observed{$rank} = $taxon;
		}
		$prevtaxon = $fulltaxon;
	}

	$prevtaxon="";
}

##########

sub reformat_lineage2 {
	if (exists $observed{"domain"}) {
		$domain = $observed{"domain"};
		$domain = $domain.";";
	}
	else {
		print "Couldn't find domain for accession $acc\n";
		$domain="unknownDomain;";
	}

	if (exists $observed{"kingdom"}) {
		$kingdom = $observed{"kingdom"};
		$kingdom = $kingdom.";";
	}
	else {
		$domain =~ s/;$//g;
		$kingdom = $domain."_undef;";
		$domain = $domain.";";
	}

	if (exists $observed{"phylum"}) {
		$phylum = $observed{"phylum"};
		$phylum = $phylum.";";

		if (exists $duplicates{$phylum}) {# handle duplicate phyla
			$kingdom =~ s/;$//;
			$phylum = $kingdom."_".$phylum;
			$kingdom = $kingdom.";";
		}
	}
	else {
		$kingdom =~ s/;$//;
		$phylum = $kingdom."_undef;";
		$kingdom = $kingdom.";";
	}

	if (exists $observed{"class"}) {
		$class = $observed{"class"};
		$class = $class.";";

		if (exists $duplicates{$class}) { #handle duplicate classes
			$phylum =~ s/;$//;
			$class = $phylum."_".$class;
			$phylum = $phylum.";";
		}
	}
	else {
		$phylum =~ s/;$//;
		$class = $phylum."_undef;";
		$phylum = $phylum.";";
	}

	if (exists $observed{"order"}) {
		$order = $observed{"order"};
		$order = $order.";";

		if (exists $duplicates{$order}) { #handle duplicate orders
			$class =~ s/;$//;
			$order = $class."_".$order;
			$class = $class.";";
		}
	}
	else {
		$class =~ s/;$//;
		$order = $class."_undef;";
		$class = $class.";";
	}

	if (exists $observed{"family"}) {
		$family = $observed{"family"};
		$family = $family.";";

		if (exists $duplicates{$family}) {#handle duplicate families
			$order =~ s/;$//;
			$family = $order."_".$family;
			$order = $order.";";
		}
	}
	else {
		$order =~ s/;$//;
		$family = $order."_undef;";
		$order = $order.";";
	}

	if (exists $observed{"genus"}) {
		$genus = $observed{"genus"};
		$genus = $genus.";";
	}

	if (exists $observed{"species"}) {
		$species = $observed{"species"};
	}

	$newlineage = $domain.$kingdom.$phylum.$class.$order.$family.$genus.$species;

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
# find the inner key of a hash of hashes with the largest value
# https://stackoverflow.com/questions/35401306/find-key-for-greatest-value-in-hash-of-hashes-in-perl
sub key_with_highest_val {
	my ($h) = @_;
	my $hi_v;
	my $hi_k;
	my $tie=0;
	for my $k (keys(%$h)) {
		my $v = $h->{$k};
		if (!defined($hi_v) || $v > $hi_v) {
			$hi_v = $v;
			$hi_k = $k;
			$tie = 0;
		}
		elsif ($v == $hi_v) {
			$tie = 1;
		}
	}
	return ($hi_k, $tie);
}

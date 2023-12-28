#!/usr/bin/perl
# Jan. 28/20 train to genus rank only, give up on resolving species issues with SILVA
#Apr. 11, 2018 retain species rank for CO1 v3git training set; handle 'superkingdom' rank and create 'cellularOrganisms'
#Jan. 5, 2017 edit to retain species rank for CO1v3 training set
#Oct. 24, 2016 add superkingdom to hash
#Aug. 11, 2016 add Protura
#March 22, 2013 by Terri Porter
#Script to create taxonomy file for Ribosomal Database Project Naive Bayesian Classifier
#usage perl make_NBC_taxonomy.plx taxid.parsed.awk.uniq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $cellularOrganisms;
my $superkingdom;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
#my $species;
my $termcounter=0;
my $j;
my $key;
my $value;
my $previous;
my $termcounter_previous;

#declare array
my @in;
my @value;
my @line;
my @previous;

#declare hash
my %lineage;
my %sorted;

open (IN, "<", $ARGV[0]) || die "Error cannot open taxid.parsed.awk.uniq: $!\n";
@in = <IN>;
close IN;

#hash lineage table and create relationships
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);
	
	$cellularOrganisms = $line[0];
	if (length $cellularOrganisms) {
		$cellularOrganisms =~ s/\s+//g; #removing trailing space
	}
	else {
		$cellularOrganisms=();
	}
	
	$superkingdom = $line[1];
	if (length $superkingdom) {
		$superkingdom =~ s/\s+//g; #remove any spaces
	}
	else {
		$superkingdom=();
	}
	
	$kingdom = $line[2];
	if (length $kingdom) {
		$kingdom =~ s/\s+//g;
	}
	else {
		$kingdom=();
	}

	$phylum = $line[3];
	if (length $phylum) {
		$phylum =~ s/\s+//g;
	}

	$class = $line[4];
	if (length $class) {
		$class =~ s/\s+//g;
	}
	else {
		$class=();
	}

	$order = $line[5];
	if (length $order) {
		$order =~ s/\s+//g;
	}
	else {
		$order=();
	}

	$family = $line[6];
	if (length $family) {
		$family =~ s/\s+//g;
	}
	else {
		$family=();
	}

	$genus = $line[7];
	if (length $genus) {
		$genus =~ s/\s+//g;
	}
	else {
		$genus=();
	}
	if (length $cellularOrganisms) {
		unless (exists $lineage{$cellularOrganisms}) {
			$termcounter++;
			$j = $termcounter-1;
			$lineage{$cellularOrganisms} = $termcounter."*".$cellularOrganisms."*".$j."*0*cellularOrganisms";
		}
	}	
	if (length $superkingdom) {
		unless (exists $lineage{$superkingdom}) {
			create_superkingdom();
		}
	}
	if (length $kingdom) {
		unless (exists $lineage{$kingdom}) {
			create_kingdom();
		}
	}
	if (length $phylum) {
		unless (exists $lineage{$phylum}) {
			create_phylum();
		}
	}
	if (length $class) {
		unless (exists $lineage{$class}) {
			create_class();
		}
	}
	if (length $order) {
		unless (exists $lineage{$order}) {
			create_order();
		}
	}		
	if (length $family) {
		unless (exists $lineage{$family}) {
			create_family();
		}
	}	
	if (length $genus) {
		unless (exists $lineage{$genus}) {
			create_genus();
		}
	}
	$i++;
	$line=();
	@line=();
	$cellularOrganisms=();
	$superkingdom=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
}
$i=0;
$termcounter=();

#sort the hash by termcounter into a new hash indexed by termcounter
while ( ($key, $value) = each (%lineage) ) {
	@value = split(/\*/, $value);
	$termcounter = $value[0];
	$sorted{$termcounter} = $value;
}
$key=();
$value=();

#sort keys then print

open (OUT, ">>", "testNBC.taxonomy") || die "Error cannot open outfile:$!\n";

foreach $key (sort {$a <=> $b} keys %sorted) {
	$value = $sorted{$key};
	print OUT $value."\n";
}
close OUT;

##########

sub create_genus {

	$termcounter++;
	if (length $family) {
		$previous = $lineage{$family};
	}
	elsif (length $order) {
		$previous = $lineage{$order};
	}
	elsif (length $class) {
		$previous = $lineage{$class};
	}
	elsif (length $phylum) {
		$previous = $lineage{$phylum};
	}
	elsif (length $kingdom) {
		$previous = $lineage{$kingdom};
	}
	elsif (length $superkingdom) {
		$previous = $lineage{$superkingdom};
	}
	else {
		$previous = $lineage{$cellularOrganisms};
	}
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$genus} = $termcounter."*".$genus."*".$termcounter_previous."*7*genus";

}

###########

sub create_family {

	$termcounter++;
	if (length $order) {
		$previous = $lineage{$order};
	}
	elsif (length $class) {
		$previous = $lineage{$class};
	}
	elsif (length $phylum) {
		$previous = $lineage{$phylum};
	}
	elsif (length $kingdom) {
		$previous = $lineage{$kingdom};
	}
	elsif (length $superkingdom) {
		$previous = $lineage{$superkingdom};
	}
	else {
		$previous = $lineage{$cellularOrganisms};
	}
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$family} = $termcounter."*".$family."*".$termcounter_previous."*6*family";

}

############

sub create_order {

	$termcounter++;
	if (length $class) {
		$previous = $lineage{$class};
	}
	elsif (length $phylum) {
		$previous = $lineage{$phylum};
	}
	elsif (length $kingdom) {
		$previous = $lineage{$kingdom};
	}
	elsif (length $superkingdom) {
		$previous = $lineage{$superkingdom};
	}
	else {
		$previous = $lineage{$cellularOrganisms};
	}	
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$order} = $termcounter."*".$order."*".$termcounter_previous."*5*order";

}

#############

sub create_class {

	$termcounter++;
	if (length $phylum) {
		$previous = $lineage{$phylum};
	}
	elsif (length $kingdom) {
		$previous = $lineage{$kingdom};
	}
	elsif (length $superkingdom) {
		$previous = $lineage{$superkingdom};
	}
	else {
		$previous = $lineage{$cellularOrganisms};
	}	
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$class} = $termcounter."*".$class."*".$termcounter_previous."*4*class";

}

##############

sub create_phylum {

	$termcounter++;
	if (length $kingdom) {
		$previous = $lineage{$kingdom};
	}
	elsif (length $superkingdom) {
		$previous = $lineage{$superkingdom};
	}
	else {
		$previous = $lineage{$cellularOrganisms};
	}	
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$phylum} = $termcounter."*".$phylum."*".$termcounter_previous."*3*phylum";

}

###############

sub create_kingdom {

	$termcounter++;
	if (length $superkingdom) {
		$previous = $lineage{$superkingdom};
	}
	else {
		$previous = $lineage{$cellularOrganisms};
	}	
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$kingdom} = $termcounter."*".$kingdom."*".$termcounter_previous."*2*kingdom";

}

###############

sub create_superkingdom {

	$termcounter++;
	if (length $cellularOrganisms) {
		$previous = $lineage{$cellularOrganisms};
	}	
	else {
		print "Can't find cellular Organisms field:\n";
	}
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$superkingdom} = $termcounter."*".$superkingdom."*".$termcounter_previous."*1*superkingdom";

}

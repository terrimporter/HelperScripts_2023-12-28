#!/usr/bin/perl
# Oct. 10/21 automate the fixing of duplicate taxon lineages
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
my $root;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
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
	
	$root = $line[0];
	$root =~ s/\s+//g; #removing trailing space
	
	$kingdom = $line[1];
	$kingdom =~ s/\s+//g;

	$phylum = $line[2];
	$phylum =~ s/\s+//g;

	$class = $line[3];
	$class =~ s/\s+//g;

	$order = $line[4];
	$order =~ s/\s+//g;

	$family = $line[5];
	$family =~ s/\s+//g;
	
	$genus = $line[6];
	$genus =~ s/\s+//g;

	$species = $line[7];
	$species =~ s/^\s+//g; #remove any preceeding whitespace
	$species =~ s/\s+/_/g; #change spaces to underscores just in case
	$species =~ s/_$//g; #if underscore at end of line, remove it

	unless (exists $lineage{$root}) {			
		$lineage{$root} = $termcounter."*".$root."*".$termcounter."*0*root";
	}

	unless (exists $lineage{$kingdom}) {
		create_kingdom();
	}

	unless (exists $lineage{$phylum}) {
		create_phylum();
	}

	unless (exists $lineage{$class}) {
		create_class();
	}

	unless (exists $lineage{$order}) {
		create_order();
	}

	unless (exists $lineage{$family}) {
		create_family();
	}

	unless (exists $lineage{$genus}) {
		create_genus();
	}

	unless (exists $lineage{$species}) {
		create_species();
	}

	$i++;
	$line=();
	@line=();
	$root=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
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

sub create_species {

	$termcounter++;
	$previous = $lineage{$genus};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$species} = $termcounter."*".$species."*".$termcounter_previous."*7*species";

}

##########

sub create_genus {

	$termcounter++;
	$previous = $lineage{$family};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$genus} = $termcounter."*".$genus."*".$termcounter_previous."*6*genus";

}

###########

sub create_family {

	$termcounter++;
	$previous = $lineage{$order};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$family} = $termcounter."*".$family."*".$termcounter_previous."*5*family";

}

############

sub create_order {

	$termcounter++;
	$previous = $lineage{$class};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$order} = $termcounter."*".$order."*".$termcounter_previous."*4*order";

}

#############

sub create_class {

	$termcounter++;
	$previous = $lineage{$phylum};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$class} = $termcounter."*".$class."*".$termcounter_previous."*3*class";

}

##############

sub create_phylum {

	$termcounter++;
	$previous = $lineage{$kingdom};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$phylum} = $termcounter."*".$phylum."*".$termcounter_previous."*2*phylum";

}

###############

sub create_kingdom {

	$termcounter++;
	$previous = $lineage{$root};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$kingdom} = $termcounter."*".$kingdom."*".$termcounter_previous."*1*kingdom";

}


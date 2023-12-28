#!/usr/bin/perl
#March 8, 2017 edit to add 'species' field for parsing CO1v3-trained to species rank
#Nov. 17, 2016 add 'phylum' field
#March 11, 2014 add 'class' field
#NEW USAGE perl count_singletons_speciestophylum.plx testNBC.fasta

#March 26, 2013 by Terri Porter
#Script to count number of orders, family, and genera represented by a single sequence in testNBC.fasta
#usage perl count_singletons.plx testNBC.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $scalar;
my $speciesField;
my $genusField;
my $familyField;
my $orderField;
my $classField;
my $phylumField;
my $species;
my $genus;
my $family;
my $order;
my $class;
my $phylum;
my $original;
my $new;
my $speciesCount=0;
my $genusCount=0;
my $familyCount=0;
my $orderCount=0;
my $classCount=0;
my $phylumCount=0;
my $count;

#declare arrays
my @in;
my @line;
my @species;
my @genus;
my @family;
my @order;
my @class;
my @phylum;

#declare hashes
my %species;
my %genus;
my %family;
my %order;
my %class;
my %phylum;

open (IN, "<", $ARGV[0]) || die "Error cannot open testNBC.fasta:$!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(";", $line);
		$scalar = scalar(@line);
		$speciesField = $scalar-1;
		$genusField = $scalar-2;
		$familyField = $scalar-3;
		$orderField = $scalar-4;
		$classField = $scalar-5;
		$phylumField = $scalar-6;
		$species = $line[$speciesField];
		$genus = $line[$genusField];
		$family = $line[$familyField];
		$order = $line[$orderField];
		$class = $line[$classField];
		$phylum = $line[$phylumField];

		#hash species, genus, family, and order names; count duplicates
		if (exists $species{$species}) {
			$original = $species{$species};
			$new = $original+1;
			$species{$species} = $new;

			$original=();
			$new=();
		}
		else {
			$species{$species} = 1;
		}

		if (exists $genus{$genus}) {
			$original = $genus{$genus};
			$new = $original+1;
			$genus{$genus} = $new;

			$original=();
			$new=();
		}
		else {
			$genus{$genus} = 1;
		}

		if (exists $family{$family}) {
			$original = $family{$family};
			$new = $original+1;
			$family{$family} = $new;

			$original=();
			$new=();
		}
		else {
			$family{$family} = 1;
		}

		if (exists $order{$order}) {
			$original = $order{$order};
			$new = $original+1;
			$order{$order} = $new;

			$original=();
			$new=();

		}
		else {
			$order{$order} = 1;
		}

		if (exists $class{$class}) {
			$original = $class{$class};
			$new = $original+1;
			$class{$class} = $new;

			$original=();
			$new=();

		}
		else {
			$class{$class} = 1;
		}

		if (exists $phylum{$phylum}) {
			$original = $phylum{$phylum};
			$new = $original+1;
			$phylum{$phylum} = $new;

			$original=();
			$new=();

		}
		else {
			$phylum{$phylum} = 1;
		}

	}
	$i++;
	$line=();
	@line=();
	$scalar=();
	$speciesField=();
	$genusField=();
	$familyField=();
	$orderField=();
	$classField=();
	$phylumField=();
	$species=();
	$genus=();
	$family=();
	$order=();
	$class=();
	$phylum=();
	
}
$i=0;

#count singletons

while ( ($species,$count) = each(%species) ) {
	if ($count == 1) {
		$speciesCount++;
		push(@species, $species);
	}
}

while ( ($genus,$count) = each(%genus) ) {
	if ($count == 1) {
		$genusCount++;
		push(@genus, $genus);
	}
}

while ( ($family,$count) = each(%family) ) {
	if ($count == 1) {
		$familyCount++;
		push(@family, $family);
	}
}

while ( ($order, $count) = each(%order) ) {
	if ($count == 1) {
		$orderCount++;
		push(@order, $order);
	}
}

while ( ($class, $count) = each(%class) ) {
	if ($count == 1) {
		$classCount++;
		push(@class, $class);
	}
}

while ( ($phylum, $count) = each(%phylum) ) {
	if ($count == 1) {
		$phylumCount++;
		push(@phylum, $phylum);
	}
}

print "\n$speciesCount singleton species\n$genusCount singleton genera\n$familyCount singleton families\n$orderCount singleton orders\n$classCount singleton classes\n$phylumCount singleton phyla\n\n";

#print singletons to files for appendix
open (OUT11, ">>", "singleton_species.txt") || die "Error cannot open singleton_species.txt: $!\n";

while ($species[$i]) {
	$species = $species[$i];
	print OUT11 "$species\n";
	$i++;
}
$i=0;
close OUT11;

open (OUT, ">>", "singleton_genera.txt") || die "Error cannot open singleton_genera.txt: $!\n";

while ($genus[$i]) {
	$genus = $genus[$i];
	print OUT "$genus\n";
	$i++;
}
$i=0;
close OUT;

open (OUT2, ">>", "singleton_families.txt") || die "Error cannot open singleton_families.txt: $!\n";

while ($family[$i]) {
	$family = $family[$i];
	print OUT2 "$family\n";
	$i++;
}
$i=0;
close OUT2;

open (OUT3, ">>", "singleton_orders.txt") || die "Error cannot open singleton_orders.txt: $!\n";

while ($order[$i]) {
	$order = $order[$i];
	print OUT3 "$order\n";
	$i++;
}
$i=0;
close OUT3;

open (OUT4, ">>", "singleton_classes.txt") || die "Error cannot open singleton_classes.txt: $!\n";

while ($class[$i]) {
	$class = $class[$i];
	print OUT4 "$class\n";
	$i++;
}
$i=0;
close OUT4;

open (OUT5, ">>", "singleton_phyla.txt") || die "Error cannot open singleton_phyla.txt: $!\n";

while ($phylum[$i]) {
	$phylum = $phylum[$i];
	print OUT5 "$phylum\n";
	$i++;
}
$i=0;
close OUT5;


#print most abundant rank categories
open (OUT12, ">>", "species_abundance.txt") || die "Error cannot open species_abundace.txt: $!\n";

foreach $species ( sort {$species{$b} <=> $species{$a}} keys %species ) {
	$count = $species{$species};
	print OUT12 "$species\t$count\n";
	$count=();
}
close OUT12;

open (OUT6, ">>", "genus_abundance.txt") || die "Error cannot open genus_abundace.txt: $!\n";

foreach $genus ( sort {$genus{$b} <=> $genus{$a}} keys %genus ) {
	$count = $genus{$genus};
	print OUT6 "$genus\t$count\n";
	$count=();
}
close OUT6;

open (OUT7, ">>", "family_abundance.txt") || die "Error cannot open family_abundance.txt: $!\n";

foreach $family ( sort {$family{$b} <=> $family{$a}} keys %family ) {
	$count = $family{$family};
	print OUT7 "$family\t$count\n";
	$count=();
}
close OUT7;

open (OUT8, ">>", "order_abundance.txt") || die "Error cannot open order_abundance.txt: $!\n";

foreach $order ( sort {$order{$b} <=> $order{$a}} keys %order ) {
	$count = $order{$order};
	print OUT8 "$order\t$count\n";
	$count=();
}
close OUT8;

open (OUT9, ">>", "class_abundance.txt") || die "Error cannot open class_abundance.txt: $!\n";

foreach $class ( sort {$class{$b} <=> $class{$a}} keys %class ) {
	$count = $class{$class};
	print OUT9 "$class\t$count\n";
	$count=();
}
close OUT9;

open (OUT10, ">>", "phylum_abundance.txt") || die "Error cannot open phylum_abundance.txt: $!\n";

foreach $phylum ( sort {$phylum{$b} <=> $phylum{$a}} keys %phylum ) {
	$count = $phylum{$phylum};
	print OUT10 "$phylum\t$count\n";
	$count=();
}
close OUT10;

#!/usr/bin/perl	
#May 13, 2013 edited to parse a testNBC.fasta that only goes to family rank
#March 26, 2013 by Terri Porter
#Script to count number of orders, family, and genera represented by a single sequence in testNBC.fasta
#usage perl count_singletons.plx testNBC.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $scalar;
my $genusField;
my $familyField;
my $orderField;
my $genus;
my $family;
my $order;
my $original;
my $new;
my $genusCount=0;
my $familyCount=0;
my $orderCount=0;
my $count;

#declare arrays
my @in;
my @line;
my @genus;
my @family;
my @order;

#declare hashes
my %genus;
my %family;
my %order;

open (IN, "<", $ARGV[0]) || die "Error cannot open testNBC.fasta:$!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(";", $line);
		$scalar = scalar(@line);
#		$genusField = $scalar-1;
		$familyField = $scalar-1;
		$orderField = $scalar-2;
#		$genus = $line[$genusField];
		$family = $line[$familyField];
		$order = $line[$orderField];

		#hash genus, family, and order names; count duplicates
#		if (exists $genus{$genus}) {
#			$original = $genus{$genus};
#			$new = $original+1;
#			$genus{$genus} = $new;

#			$original=();
#			$new=();
#		}
#		else {
#			$genus{$genus} = 1;
#		}

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
	}
	$i++;
	$line=();
	@line=();
	$scalar=();
#	$genusField=();
	$familyField=();
	$orderField=();
#	$genus=();
	$family=();
	$order=();
	
}
$i=0;

#count singletons
#while ( ($genus,$count) = each(%genus) ) {
#	if ($count == 1) {
#		$genusCount++;
#		push(@genus, $genus);
#	}
#}

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

print "\n$familyCount singleton families\n$orderCount singleton orders\n\n";

#print singletons to files for appendix
#open (OUT, ">>", "singleton_genera.txt") || die "Error cannot open singleton_genera.txt: $!\n";

#while ($genus[$i]) {
#	$genus = $genus[$i];
#	print OUT "$genus\n";
#	$i++;
#}
#$i=0;
#close OUT;

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

#print most abundant rank categories

#open (OUT4, ">>", "genus_abundance.txt") || die "Error cannot open genus_abundace.txt: $!\n";

#foreach $genus ( sort {$genus{$b} <=> $genus{$a}} keys %genus ) {
#	$count = $genus{$genus};
#	print OUT4 "$genus\t$count\n";
#	$count=();
#}
#close OUT4;

open (OUT5, ">>", "family_abundance.txt") || die "Error cannot open family_abundance.txt: $!\n";

foreach $family ( sort {$family{$b} <=> $family{$a}} keys %family ) {
	$count = $family{$family};
	print OUT5 "$family\t$count\n";
	$count=();
}
close OUT5;

open (OUT6, ">>", "order_abundance.txt") || die "Error cannot open order_abundance.txt: $!\n";

foreach $order ( sort {$order{$b} <=> $order{$a}} keys %order ) {
	$count = $order{$order};
	print OUT6 "$order\t$count\n";
	$count=();
}
close OUT6;

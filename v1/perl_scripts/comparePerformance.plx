#!/usr/bin/perl
#March 21, 2012 by Terri Porter
#Script to compare results using RDP with different reference libraries trained for different ranks
#usage perl comparePerformance.plx features.txt.parsed_old output_rank_LR3

use strict;
use warnings;

#declare var
my $RANK="family"; #####specify rank (lowercase) to provide performance stats for (genus - kingdom)#####
my $CUTOFF=0.5; #####specify cutoff (0-1) here#####
my $i=0;
my $line;
my $gb;
my $binomial;
my $lineage;
my $j=0;
my $gb_test;
my $kingdom_test;
my $kingdom_test_boot;
my $phylum_test;
my $phylum_test_boot;
my $class_test;
my $class_test_boot;
my $order_test;
my $order_test_boot;
my $family_test;
my $family_test_boot;
my $genus_test;
my $genus_test_boot;
my $boot;
my $ref;
my $test;
my $recovery=0;
my $coverage=0;
my $error=0;
my $count=0;
my $key;
my $value;
my $genus;
my $name;

#declare array
my @features;
my @line;
my @lineage;
my @rdp;
my @binomial;

#declare hash
my %kingdom;
my %phylum;
my %class;
my %order;
my %family;
my %genus;
my %kingdom_test;
my %phylum_test;
my %class_test;
my %order_test;
my %family_test;
my %genus_test;
my %kingdom_test_boot;
my %phylum_test_boot;
my %class_test_boot;
my %order_test_boot;
my %family_test_boot;
my %genus_test_boot;

open (FEATURES,"<",$ARGV[0]) || die "Error reading features.txt.parsed_old: $!\n";
@features = <FEATURES>;
close FEATURES;

#populate rank hashes with reference names for correct assignments indexed by gb accession

while ($features[$i]) {
	$line = $features[$i];
	chomp $line;

	if ($line =~ /^\S+/) {
		@line = split(/\t/,$line);
		$gb = $line[0];
		$binomial = $line[1];
		@binomial = split(/\s/,$binomial);
		$genus = $binomial[0];
		
		$lineage = $line[5];
		@lineage = split(/\s/,$lineage);
#		$genus = pop(@lineage);

		while($lineage[$j]) {
			$name = $lineage[$j];
			
			if ($name =~ /Fungi/) {
				$kingdom{$gb} = "Fungi";
			}
			elsif ($name =~ /mycota$/) {
				$phylum{$gb} = $name;
			}
			elsif ($name =~ /etes$/) {
				$class{$gb} = $name;
			}
			elsif ($name =~ /ales$/) {
				$order{$gb} = $name;
			}
			elsif ($name =~ /ceae$/) {
				$family{$gb} = $name;
			}
			elsif ($name eq $genus) {
				$genus{$gb} = $name;
			}
			$j++;
			$name=();
		}
		$j=0;
	}
	$i++;
	@line=();
	$gb=();
#	$binomial=();
	$lineage=();
	@lineage=();
	$genus=();

}
$i=0;

#test
while(($key,$value) = each(%genus)) {
	print "genus key->$key\t value->$value\n";
}

open (RDP,"<",$ARGV[1]) || die "Error cannot read rdp_output.txt: $!\n";
@rdp = <RDP>;
close RDP;

#populate rank hases for comparison

while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	if ($line =~ /^\S+/	) {
		@line = split(/\t/,$line);
		$gb_test = $line[0];
		$kingdom_test = $line[5];
		$kingdom_test{$gb_test} = $kingdom_test;
		$kingdom_test_boot = $line[7];
		$kingdom_test_boot{$gb_test} = $kingdom_test_boot;
		$phylum_test = $line[8];
		$phylum_test{$gb_test} = $phylum_test;
		$phylum_test_boot = $line[10];
		$phylum_test_boot{$gb_test} = $phylum_test_boot;
		$class_test = $line[11];
		$class_test{$gb_test} = $class_test;
		$class_test_boot = $line[13];
		$class_test_boot{$gb_test} = $class_test_boot;
		$order_test = $line[14];
		$order_test{$gb_test} = $order_test;
		$order_test_boot = $line[16];
		$order_test_boot{$gb_test} = $order_test_boot;
		$family_test = $line[17];
		$family_test{$gb_test} = $family_test;
		$family_test_boot = $line[19];
		$family_test_boot{$gb_test} = $family_test_boot;
		$genus_test = $line[20];
		$genus_test{$gb_test} = $genus_test;
		$genus_test_boot = $line[22];
		$genus_test_boot{$gb_test} = $genus_test_boot;
	}
	$i++;
	@line=();
	$gb_test=();
	$kingdom_test=();
	$kingdom_test_boot=();
	$phylum_test=();
	$phylum_test_boot=();
	$class_test=();
	$class_test_boot=();
	$order_test=();
	$order_test_boot=();
	$family_test=();
	$family_test_boot=();
	$genus_test=();
	$genus_test_boot=();
}
$i=0;

#compare recovery at specified rank with specified cutoff

if ($RANK eq "genus") {

	while( ($key,$value)  = each(%genus_test_boot) ) {
		$gb_test = $key;
		$boot = $value;
		
		if ($boot >= $CUTOFF) {

			$ref = $genus{$gb_test};
			$test = $genus_test{$gb_test};

			if ($ref eq $test) {
				$recovery++;
				$coverage++;
			}
			else {
				$error++;
				$coverage++;
			}
		}
		$count++;
	}
	print "Recovery: $recovery\nErroroneous recovery: $error\nCoverage: $coverage\nTotal count: $count\n";
	$recovery=0;
	$error=0;
	$coverage=0;
	$count=0;

}

elsif ($RANK eq "family") {

	while( ($key,$value)  = each(%family_test_boot)) {
		$gb_test = $key;
		$boot = $value;
		
		if ($boot >= $CUTOFF) {

			$ref = $family{$gb_test};
			$test = $family_test{$gb_test};

			if ($ref eq $test) {
				$recovery++;
				$coverage++;
			}
			else {
				$error++;
				$coverage++;
			}
		}
		$count++;
	}
	print "Recovery: $recovery\nErroroneous recovery: $error\nCoverage: $coverage\nTotal count: $count\n";
	$recovery=0;
	$error=0;
	$coverage=0;
	$count=0;
}

elsif ($RANK eq "order") {

	while( ($key,$value)  = each(%order_test_boot)) {
		$gb_test = $key;
		$boot = $value;
		
		if ($boot >= $CUTOFF) {

			$ref = $order{$gb_test};
			$test = $order_test{$gb_test};

			if ($ref eq $test) {
				$recovery++;
				$coverage++;
			}
			else {
				$error++;
				$coverage++;
			}
		}
		$count++;
	}
	print "Recovery: $recovery\nErroroneous recovery: $error\nCoverage: $coverage\nTotal count: $count\n";
	$recovery=0;
	$error=0;
	$coverage=0;
	$count=0;
}

elsif ($RANK eq "class") {

	while( ($key,$value)  = each(%class_test_boot)) {
		$gb_test = $key;
		$boot = $value;
		
		if ($boot >= $CUTOFF) {

			$ref = $class{$gb_test};
			$test = $class_test{$gb_test};

			if ($ref eq $test) {
				$recovery++;
				$coverage++;
			}
			else {
				$error++;
				$coverage++;
			}
		}
		$count++;
	}
	print "Recovery: $recovery\nErroroneous recovery: $error\nCoverage: $coverage\nTotal count: $count\n";
	$recovery=0;
	$error=0;
	$coverage=0;
	$count=0;
}

elsif ($RANK eq "phylum") {

	while( ($key,$value)  = each(%phylum_test_boot)) {
		$gb_test = $key;
		$boot = $value;
		
		if ($boot >= $CUTOFF) {

			$ref = $phylum{$gb_test};
			$test = $phylum_test{$gb_test};

			if ($ref eq $test) {
				$recovery++;
				$coverage++;
			}
			else {
				$error++;
				$coverage++;
			}
		}
		$count++;
	}
	print "Recovery: $recovery\nErroroneous recovery: $error\nCoverage: $coverage\nTotal count: $count\n";
	$recovery=0;
	$error=0;
	$coverage=0;
	$count=0;
}

elsif ($RANK eq "kingdom") {

	while( ($key,$value)  = each(%kingdom_test_boot)) {
		$gb_test = $key;
		$boot = $value;
		
		if ($boot >= $CUTOFF) {

			$ref = $kingdom{$gb_test};
			$test = $kingdom_test{$gb_test};

			if ($ref eq $test) {
				$recovery++;
				$coverage++;
			}
			else {
				$error++;
				$coverage++;
			}
		}
		$count++;
	}
	print "Recovery: $recovery\nErroroneous recovery: $error\nCoverage: $coverage\nTotal count: $count\n";
	$recovery=0;
	$error=0;
	$coverage=0;
	$count=0;
}


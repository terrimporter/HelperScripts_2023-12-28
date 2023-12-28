#!/usr/bin/perl
#July 8, 2013 edit to work with Mantodea dataset with v4 GenBank-family trained classifier
#March 28, 2013 edit to grab stats for all fields (class to genus)
#March 25, 2013 by Terri Porter
#Script to check classification tests using subusets of testNBC.fasta classified using testNBC.taxonomy
#usage perl check_NBC.plx testquery.out testNBC.fasta
#add a filter to only look at assignments with > x% bootstrap confidence

use strict;
use warnings;

#declare vars
my $i=0;
my $line;
my $scalar;
my $genusField;
my $id;
my $genusAsst;
my $idstring;
my $bootstrapCutoff = 0.96; ### Edit minimum bootstrap cutoff here if needed ###
my $genusBootstrapField;
my $genusConfidence;
my $familyField;
my $familyBootstrapField;
my $familyAsst;
my $familyConfidence;
my $orderField;
my $orderBootstrapField;
my $orderAsst;
my $orderConfidence;
my $classField;
my $classBootstrapField;
my $classAsst;
my $classConfidence;
my $asstString;
my $confidenceString;
my $knownString;
my $genusCorrect=0;
my $genusWrong=0;
my $genusNotClassified=0;
my $familyCorrect=0;
my $familyWrong=0;
my $familyNotClassified=0;
my $orderCorrect=0;
my $orderWrong=0;
my $orderNotClassified=0;
my $classCorrect=0;
my $classWrong=0;
my $classNotClassified=0;
my $classKnown;
my $orderKnown;
my $familyKnown;
my $genusKnown;
my $classCountsRef;

#declare arrays
my @nbc;
my @line;
my @known;
my @asstString;
my @confidenceString;
my @knownString;
my @classCounts = (0,0,0); #(correct, wrong, not classified)
my @orderCounts = (0,0,0);
my @familyCounts = (0,0,0);
my @genusCounts = (0,0,0);

#print "genusCounts: @genusCounts\n";

#declare hashes
my %asst;
my %known;
my %confidence;

open (NBC, "<", $ARGV[0]) || die "Error cannot open testquery.out: $!\n";
@nbc = <NBC>;
close NBC;

open (KNOWN, "<", $ARGV[1]) || die "Error cannot open testNBC.fasta: $!\n";
@known = <KNOWN>;
close KNOWN;

#hash assts
while ($nbc[$i]) {
	$line = $nbc[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$scalar = scalar(@line);
	$id = $line[0];

	#genus
#	$genusField = $scalar-3;
#	$genusBootstrapField = $scalar-1;
#	$genusAsst = $line[$genusField];
#	$genusConfidence = $line[$genusBootstrapField];

	#family
	$familyField = $scalar-3;
	$familyBootstrapField = $scalar-1;
	$familyAsst = $line[$familyField];
	$familyConfidence = $line[$familyBootstrapField];

	#order
	$orderField = $scalar-6;
	$orderBootstrapField = $scalar-4;
	$orderAsst = $line[$orderField];
	$orderConfidence = $line[$orderBootstrapField];

	#class Insecta
	$classField = $scalar-9;
	$classBootstrapField = $scalar-7;
	$classAsst = $line[$classField];
	$classConfidence = $line[$classBootstrapField];

#	$asst{$id} = $classAsst.";".$orderAsst.";".$familyAsst.";".$genusAsst;
	$asst{$id} = $classAsst.";".$orderAsst.";".$familyAsst;
#	$confidence{$id} = $classConfidence.";".$orderConfidence.";".$familyConfidence.";".$genusConfidence;
	$confidence{$id} = $classConfidence.";".$orderConfidence.";".$familyConfidence;

	$i++;
	$line=();
	@line=();
	$scalar=();
	$genusField=();
	$genusBootstrapField=();
	$genusAsst=();
	$genusConfidence=();
	$familyField=();
	$familyBootstrapField=();
	$familyAsst=();
	$familyConfidence=();
	$orderField=();
	$orderBootstrapField=();
	$orderAsst=();
	$orderConfidence=();
	$classField=();
	$classBootstrapField=();
	$classAsst=();
	$classConfidence=();
	$id=();
}
$i=0;

#hash known
while ($known[$i]) {
	$line = $known[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(";", $line);
		$idstring = $line[0];
		$idstring =~ /^>(\S+)\s+/;
		$id = $1;
		$scalar = scalar(@line);
		$genusField = $scalar-1;
		$genusAsst = $line[$genusField];
		$familyField = $scalar-2;
		$familyAsst = $line[$familyField];
		$orderField = $scalar-3;
		$orderAsst = $line[$orderField];
		$classField = $scalar-4;
		$classAsst = $line[$classField];
		
#		$known{$id} = $classAsst.";".$orderAsst.";".$familyAsst.";".$genusAsst;
		$known{$id} = $classAsst.";".$orderAsst.";".$familyAsst;

	}
	$i++;
	$line=();
	@line=();
	$idstring=();
	$id=();
	$scalar=();
	$genusField=();
	$genusAsst=();
	$familyField=();
	$familyAsst=();
	$orderField=();
	$orderAsst=();
	$classField=();
	$classAsst=();
}
$i=0;

#compare asst with known
while ( ($id, $asstString) = each %asst) {
	@asstString = split(/;/, $asstString);

	if (exists $confidence{$id}) {
		$confidenceString = $confidence{$id};
		@confidenceString = split(/;/, $confidenceString);

		if (exists $known{$id}) {
			$knownString = $known{$id};
			@knownString = split(/;/, $knownString);

			#process class
			$classAsst = $asstString[0];
			$classConfidence = $confidenceString[0];
			$classKnown = $knownString[0];
			@classCounts = compare(\$classAsst,\$classConfidence,\$classKnown,\@classCounts);

			#process order
			$orderAsst = $asstString[1];
			$orderConfidence = $confidenceString[1];
			$orderKnown = $knownString[1];
			@orderCounts = compare(\$orderAsst,\$orderConfidence,\$orderKnown,\@orderCounts);

			#process family
			$familyAsst = $asstString[2];
			$familyConfidence = $confidenceString[2];
			$familyKnown = $knownString[2];
			@familyCounts = compare(\$familyAsst,\$familyConfidence,\$familyKnown,\@familyCounts);

			#process genus
#			$genusAsst = $asstString[3];
#			$genusConfidence = $confidenceString[3];
#			$genusKnown = $knownString[3];
#			@genusCounts = compare(\$genusAsst,\$genusConfidence,\$genusKnown,\@genusCounts);

		}
		else {
			print "Error missing a known classification field for $id\n";
		}
	}
	else {
		print "Error missing bootstrap field for $id\n";
	}

	@asstString=();
	$confidenceString=();
	@confidenceString=();
	@knownString=();
	$classAsst=();
	$classConfidence=();
	$classKnown=();
	$orderAsst=();
	$orderConfidence=();
	$orderKnown=();
	$familyAsst=();
	$familyConfidence=();
	$familyKnown=();
	$genusAsst=();
	$genusConfidence=();
	$genusKnown=();

}

#print "0.9 Bootstrap cutoff\tClass\tOrder\tFamly\tGenus\n";
print "0.9 Bootstrap cutoff\tClass\tOrder\tFamly\n";

#print "Correct\t$classCounts[0]\t$orderCounts[0]\t$familyCounts[0]\t$genusCounts[0]\n";
print "Correct\t$classCounts[0]\t$orderCounts[0]\t$familyCounts[0]\n";
#print "Wrong\t$classCounts[1]\t$orderCounts[1]\t$familyCounts[1]\t$genusCounts[1]\n";
print "Wrong\t$classCounts[1]\t$orderCounts[1]\t$familyCounts[1]\n";
#print "Not classified\t$classCounts[2]\t$orderCounts[2]\t$familyCounts[2]\t$genusCounts[2]\n";
print "Not classified\t$classCounts[2]\t$orderCounts[2]\t$familyCounts[2]\n";

####################

sub compare {

my ($asstRef,$confidenceRef,$knownRef,$arrayRef) = @_;
my $asst = ${$asstRef};
my $confidence = ${$confidenceRef};
my $known = ${$knownRef};
my @counts = @{$arrayRef};

my $correct = $counts[0];
my $wrong = $counts[1];
my $notClassified = $counts[2];
my $correctNEW;
my $wrongNEW;
my $notClassifiedNEW;

	if ($confidence > $bootstrapCutoff) {
		if ($asst eq $known) {
			$correctNEW = $correct+1;
			$counts[0] = $correctNEW;
		}
		else {
			$wrongNEW = $wrong+1;
			$counts[1] = $wrongNEW;
		}
	}
	else {
		$notClassifiedNEW = $notClassified+1;
		$counts[2] = $notClassifiedNEW;
	}

	return (@counts);

}

####################

#!/usr/bin/perl
#Terri Porter by Feb. 25, 2013
#Forward simulation of tree ((A,B),C), t0 = initial population, t1 = first split, t2 = one lineage splits again, t3 = current sample
#usage perl forward.plx

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare global variables
my $resample=100; #number of simulations to run
my $theta=4;
my $popsize = 100;
my $t0t1 = $popsize/10; #divergence times
my $t1t2 = $popsize/10;
my $t2t3 = $popsize/10;
my $samplesizeA=0;#sample sizes from each population
my $samplesizeB=1;
my $samplesizeC=1;
my $length = 100;
my $mu = ($theta/(4*$popsize)); #mutation rate per generation (diploid)

my $filename;
my $i=0;
my $base;
my $sequence;
my $pattern;
my $count;
my $proportion;

my $segregatingsitesAA;
my $segregatingsitesBB;
my $segregatingsitesCC;
my $segregatingsitesAB;
my $segregatingsitesAC;
my $segregatingsitesBC;
my $samplesizeAB=$samplesizeA+$samplesizeB;
my $samplesizeAC=$samplesizeA+$samplesizeC;
my $samplesizeBC=$samplesizeB+$samplesizeC;
my $segregatingsites;

#declare global arrays
my @population_t0; #population sequences identical
my @segregatingsites; #from one simulation
my @population_t3_A_sample; #subsample of population at t3
my @population_t3_B_sample;
my @population_t3_C_sample;
my @population_t3_AB_sample;
my @population_t3_AC_sample;
my @population_t3_BC_sample;
my @allsegregatingsites; #array of @segregatingsites arrays
my @segregatingsitesElement;

#declare global hash
my %pattern;

#create custom outfile
$filename = $samplesizeA."_".$samplesizeB."_".$samplesizeC.".txt";
open (OUT, ">>",$filename) || die "Error cannot open outfile: $!\n";

#print initial values
print OUT "Theta\tPopSize\tt0t1\tt1t2\tt2t3\tgenelength\tsamplesizeA\tsamplesizeB\tsamplesizeC\tresample\n";
print OUT "$theta\t$popsize\t$t0t1\t$t1t2\t$t2t3\t$length\t$samplesizeA\t$samplesizeB\t$samplesizeC\t$resample\n\n";
print OUT "within populations:\t\t\tamong populations:\t\t\n";
print OUT "AA\tBB\tCC\tAB\tAC\tBC\n";

#create a sequence to clone
while ($i<$length) {
	$base = randomNucleotide();
	$sequence = $sequence.$base;
	$i++;
}
#print "ORIGINAL SEQUENCE\n$sequence\n\n";
$i=0;

#create initial population of identical sequences
while ($i<$popsize) {
	push(@population_t0,$sequence);
	$i++;
}
$i=0;
#print "POPULATION t0:\n";
#foreach my $seq (@population_t0) {
#	print $seq."\n";
#}
#print "\n\n";

#run simulation multiple times and add patterns to hash
while ($i<$resample) {

	@segregatingsites = simulation(\@population_t0); #[0]AA,BB,CC,AB,AC,[5]BC
#	$allsegregatingsites[$i] = [@segregatingsites]; #build array of arrays indexed by resamplenumber-1
	$pattern = join('_',@segregatingsites);
	if (exists $pattern{$pattern}) {
		$count = $pattern{$pattern};
		$count+=1;
		$pattern{$pattern} = $count;
	}
	else {
		$pattern{$pattern} = 1;
	}

	@segregatingsites=();
	$i++;
}
$i=0;

#print proportion for segregating site patterns to identify the most abundant, print most abundant firsti
print OUT "segregating site pattern\tfrequency\tproportion\n";
foreach $pattern (sort { $pattern{$b} <=> $pattern{$a}} (keys %pattern)) {
	$count = $pattern{$pattern};
	$proportion = $count/$resample;
	print OUT "$pattern\t$count\t$proportion\n";
	$count=();
	$pattern=();
	$proportion=();
}
close OUT;

##########

sub simulation {

my $j=0;
my @population_t0 = @{$_[0]};
my @population_t1_ABC; #population sequences at t1
my @population_t1_AB;
my @population_t1_C;
my @population_t2_AB;
my @population_t2_A; #population sequences at t2
my @population_t2_B;
my @population_t2_C;
my @population_t3_A; #population sequences at t3
my @population_t3_B;
my @population_t3_C;
my @population_t3_AB;
my @population_t3_AC;
my @population_t3_BC;

#t0 population evolves to t1
while ($j<$t0t1) { 
	@population_t0 = randomlyMate(\@population_t0);
	@population_t1_ABC = evolvePopulation(\@population_t0);
	$j++;
}
$j=0;

@population_t1_AB = randomlyMate(\@population_t1_ABC);
@population_t1_C = randomlyMate(\@population_t1_ABC);

#print "POPULATION t1:\n";
#foreach my $seq (@population_t1_AB) {
#	print $seq."\n";
#}
#print "\n\n";

#t1_AB population evolves to t2
while ($j<$t1t2) {
	@population_t2_AB = evolvePopulation(\@population_t1_AB);
	$j++;
}
$j=0;

@population_t2_A = randomlyMate(\@population_t2_AB);
@population_t2_B = randomlyMate(\@population_t2_AB);

#print "POPULATION t2_AB:\n";
#foreach my $seq (@population_t2_A) {
#	print $seq."\n";
#}
#print "\n\n";

#t1_C population evolves t2
while ($j<$t1t2) {
	@population_t2_C = evolvePopulation(\@population_t1_C);
	$j++;
}
$j=0;

#print "POPULATION t2_C:\n";
#foreach my $seq (@population_t2_C) {
#	print $seq."\n";
#}
#print "\n\n";

#t2_A population evolves to t3_A
while ($j<$t2t3) {
	@population_t3_A = evolvePopulation(\@population_t2_A);
	$j++;
}
$j=0;

#print "POPULATION t3_A:\n";
#foreach my $seq (@population_t3_A) {
#	print $seq."\n";
#}
#print "\n\n";

#t2_B population evolves to t3_B
while ($j<$t2t3) {
	@population_t3_B = evolvePopulation(\@population_t2_B);
	$j++;
}
$j=0;

#print "POPULATION t3_B:\n";
#foreach my $seq (@population_t3_B) {
#	print $seq."\n";
#}
#print "\n\n";

#t2_C population evolves to t3_C
while ($j<$t2t3) {
	@population_t3_C = evolvePopulation(\@population_t2_C);
	$j++;
}
$j=0;

#print "POPULATION t3_C:\n";
#foreach my $seq (@population_t3_C) {
#	print $seq."\n";
#}
#print "\n\n";

#at t3 subsample population, count number of segregating sites
@population_t3_A_sample = subsamplePopulation(\@population_t3_A,\$samplesizeA);
$segregatingsitesAA = countSegregatingSites(\@population_t3_A_sample);
push(@segregatingsites,$segregatingsitesAA);

@population_t3_B_sample = subsamplePopulation(\@population_t3_B,\$samplesizeB);
$segregatingsitesBB = countSegregatingSites(\@population_t3_B_sample);
push(@segregatingsites, $segregatingsitesBB);

@population_t3_C_sample = subsamplePopulation(\@population_t3_C,\$samplesizeC);
$segregatingsitesCC = countSegregatingSites(\@population_t3_C_sample);
push(@segregatingsites, $segregatingsitesCC);

if ($samplesizeA > 0 && $samplesizeB > 0) {
	push(@population_t3_AB_sample, @population_t3_A_sample);
	push(@population_t3_AB_sample, @population_t3_B_sample);
	$segregatingsitesAB = countSegregatingSites(\@population_t3_AB_sample);
	push(@segregatingsites, $segregatingsitesAB);
}
else {
	$segregatingsitesAB=0;
	push(@segregatingsites, $segregatingsitesAB);
}

if ($samplesizeA > 0 && $samplesizeC > 0) {
	push(@population_t3_AC_sample, @population_t3_A_sample);
	push(@population_t3_AC_sample, @population_t3_C_sample);
	$segregatingsitesAC = countSegregatingSites(\@population_t3_AC_sample);
	push(@segregatingsites, $segregatingsitesAC);
}
else {
	$segregatingsitesAC = 0;
	push (@segregatingsites, $segregatingsitesAC);
}

if ($samplesizeB > 0 && $samplesizeC > 0) {
	push(@population_t3_BC_sample, @population_t3_B_sample);
	push(@population_t3_BC_sample, @population_t3_C_sample);
	$segregatingsitesBC = countSegregatingSites(\@population_t3_BC_sample);
	push(@segregatingsites, $segregatingsitesBC);
}
else {
	$segregatingsitesBC = 0;
	push(@segregatingsites, $segregatingsitesBC);
}

print OUT "$segregatingsitesAA\t$segregatingsitesBB\t$segregatingsitesCC\t$segregatingsitesAB\t$segregatingsitesAC\t$segregatingsitesBC\n";

@population_t0=();
@population_t1_AB=();
@population_t1_C=();
@population_t2_A=();
@population_t2_B=();
@population_t2_C=();
@population_t3_A=();
@population_t3_B=();
@population_t3_C=();
@population_t3_A_sample=();
@population_t3_B_sample=();
@population_t3_C_sample=();
@population_t3_AB_sample=();
@population_t3_AC_sample=();
@population_t3_BC_sample=();
$segregatingsitesAA=();
$segregatingsitesBB=();
$segregatingsitesCC=();
$segregatingsitesAB=();
$segregatingsitesAC=();
$segregatingsitesBC=();

return @segregatingsites;

}

##########

sub randomlyMate { #randomly sample WITH replacement

my @original_population = @{$_[0]};
my @new_population=();
my $o=0;
my $individual;

while ($o < $popsize) {
	$individual = $original_population[rand @original_population];
	push(@new_population, $individual);
	$o++;
	$individual=();
}
$o=0;

@original_population=();
return @new_population;

}

##########

sub subsamplePopulation {

my @original_population = @{$_[0]};
my $sample;
my @population_sample=();
my $samplesize = ${$_[1]};

my $n=0;

while ($n < $samplesize) {
	$sample = $original_population[rand @original_population];
	push(@population_sample, $sample);
	$n++;
	$sample=();
}
$n=0;
@original_population=();
$samplesize=();
return @population_sample;

}

##########

sub countSegregatingSites {

my $k=0;
my $l=0;
my $alignmentline;
my $alignmentbase;
my $m=0;
my $original;
my $segregatingsite=0; #flag
my $sscounter=0;

my @alignment = @{$_[0]};
my @column;

while ($k< $length) {
	while ($alignment[$l]) {
		$alignmentline = $alignment[$l];
		$alignmentbase = substr($alignmentline,$k,1);
		push(@column,$alignmentbase);
		$l++;
	}
	$l=0;

	while ($column[$m]) {
		if ($m==0) {
			$original = $column[$m];
		}
		elsif ($m>0) {
			$alignmentbase = $column[$m];
			if ($original eq $alignmentbase) {
				$m++;
				next;
			}
			elsif ($original ne $alignmentbase) {
				$segregatingsite = 1;
			}
		}
		$m++;
		$alignmentbase=();
	}
	$m=0;
	@column=();
	$original=();

	if ($segregatingsite==1) {
		$sscounter++;
	}
	$segregatingsite=0;
	$k++;
	$alignmentline=();
	$alignmentbase=();
}
$k=0;

return $sscounter;

}

##########

sub evolvePopulation{

my $sequence;
my @sequence;
my $randomNumber;
#my $position;
my $base;
my $newbase;
my @mutating = @{$_[0]};
my $k=0;
my $l=0;

	while($mutating[$k]) { #flip through the sequences
		$sequence = $mutating[$k];
#		print "sequence: $sequence\n";#test
		@sequence = split(//,$sequence);
#		print "@ sequence: @sequence\n";

		while ($sequence[$l]) { #flip through the bases
			$base = $sequence[$l];
			$randomNumber = rand(); #random between 0 and 1 ############ HERE

			if ($randomNumber <= $mu) {

				if ($base eq "A") {
					$newbase = changeA();
				}
				elsif ($base eq "C") {
					$newbase = changeC();
				}
				elsif ($base eq "G") {
					$newbase = changeG();
				}
				elsif ($base eq "T") {
					$newbase = changeT();
				}
				$sequence[$l] = $newbase;
			}
			$l++;
			$randomNumber=();
			$base=();
			$newbase=();
		}
		$l=0;
		$sequence = join('',@sequence);
#		print "new sequence: $sequence\n\n";
		$mutating[$k] = $sequence;
		$k++;
		$sequence=();
		@sequence=();
	}
	$k=0;

	return @mutating;

}

##########

sub randomPosition {

my $sequence = ${$_[0]};
	return int(rand(length $sequence));

}

##########

sub randomNucleotide {

#my @DNA = ("A","C","G","T");
my @DNA = ("A","A","A","A");### testing ###

	return $DNA[rand @DNA];

}

##########

sub changeA {

my @changeA = ("C","G","T");
	
	return $changeA[rand @changeA];

}

##########

sub changeC {

my @changeC = ("A","G","T");

	return $changeC[rand @changeC];

}

##########

sub changeG {

my @changeG = ("A","C","T");

	return $changeG[rand @changeG];

}

##########

sub changeT {

my @changeT = ("A","C","G");

	return $changeT[rand @changeT];

}

##########

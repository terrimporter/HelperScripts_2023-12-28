#!/usr/bin/perl
# Teresita M. Porter, Feb. 4, 2020
# Script to introduce mutations into BOLD seqs for simulation study on pseudogenes
# usage perl mutate.plx bold.nt.fasta

use strict;
use warnings;
use Data::Dumper;
use Math::Round;

# params to adjust
my $targetSeqs = 100000; # target community size; used with seqCounter
my $mutateProp = 0.19; # proportion of library sequence to mutate; use with numtCounter and targetNumts
my $reduceGC = 0.025; # reduce GC content (%); use with counter and target
my $indelProp = 0.025; # proportion of indels to introduce; use with counter and target
my $outfile0 = "mutated_0.fasta";
my $outfile1 = "mutated_1.fasta";
my $outfile2 = "mutated_2.fasta";

# vars
my $i=0;
my $line;
my $header;
my $seq;
my $acc;
my $length;
my $x; # random location on seq
my $y; # random choice from substitution array
my $z; # random choice from substitution array2
my $base;
my $newbase;
my $counter=0;
my $target;
my $hashsize;
my $targetNumts;
my $numtCounter=0;
my $seqCounter=0;

# array
my @in;
my @line;
my @seq;
my @AT = ("A","T"); # substitution array
my @ACGT = ("A", "C", "G", "T");

# hashes
my %seqs; # key = accession, value = seq
my %sim0; # key = accession, value = seq
my %sim1; # key = accession, value = seq
my %sim2; # key = accession, value = seq

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@in = <IN>;
close IN;

# hash seqs
# expecting strict FASTA format
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$header = $line;
		@line = split(/ /, $line);
		$acc = shift @line;
	}
	else {
		$seq = $line;
		$seqs{$acc} = $seq;

	}
	$i++;
}
$i=0;

##################################################################
# simulation 0
# create a mock community (no mutations)

print "Starting to prepare outfile0\n";
open (OUT0, ">>", $outfile0) || die "Cannot open outfile0: $!\n";

while ( ($acc, $seq) = each (%seqs) ) {

	if ($seqCounter < $targetSeqs) {
		$seq = $seqs{$acc};
		$sim0{$acc} = $seq;
		$seqCounter++;
		print OUT0 ">".$acc."\n".$seq."\n";
	}
	else {
		last;
	}
}

close OUT0;

##################################################################
# simulation 1
# reduce GC composition
# use community sampled in simulation 0

# calculate proportion of library sequences to mutate
$targetNumts = int($targetSeqs * $mutateProp);


while ( ($acc, $seq) = each (%sim0) ) {

	if ($numtCounter < $targetNumts) {
		$seq = $sim0{$acc};
		@seq = split(//,$seq);
		$length = scalar @seq;
		# target is number of bases to mutate 
		$target = int($length * $reduceGC);

		reduceGC();

		$acc = $acc."_pseudogene";
		$sim1{$acc} = $seq; # replace with the mutated sequence
		$numtCounter++;
	}
	else {
		$seq = $sim0{$acc};
		$sim1{$acc} = $seq;
	}
}

print "Starting to prepare outfile1\n";

open (OUT1, ">>", $outfile1) || die "Cannot open outfile1: $!\n";

while ( ($acc, $seq)  = each (%sim1) ) {
	$seq = $sim1{$acc};
	print OUT1 ">".$acc."\n".$seq."\n";
}

close OUT1;

# reset counter
$numtCounter=0;

##################################################################
# start mutation strategy 2
# introduce indels
# use community sampled in simulation 0

while ( ($acc, $seq) = each (%sim0) ) {

	if ($numtCounter < $targetNumts) {
		$seq = $sim0{$acc};
		@seq = split(//,$seq);
		$length = scalar @seq;	
		$target = int($length * $indelProp);
			
		addIndels();

		$acc = $acc."_pseudogene";
		$sim2{$acc} = $seq; # replace with the mutated sequence
		$numtCounter++;
	}
	else {
		$seq = $seqs{$acc};
		$sim2{$acc} = $seq;
	}
}

print "Starting to prepare outfile2\n";

open (OUT2, ">>", $outfile2) || die "Cannot open outfile2: $!\n";

while ( ($acc, $seq)  = each (%sim2) ) {
	$seq = $sim2{$acc};
	print OUT2 ">".$acc."\n".$seq."\n";
}

close OUT2;

###############################

sub reduceGC {

	# target is number of bases to mutate
	while ($counter <= $target ) {

		# pick a random position along sequence 
		# search between 0 - length (doesn't include length val)
		$x = int(rand($length));
		$base = substr($seq, $x, 1);

		# replace one G/C with an A/T
		if ($base =~ /[CG]/) {
			# rand picks fractional number between 0 < 1 
			# rounds number to nearest integer
			$y = round(rand(1));
			$newbase = $AT[$y];
			substr($seq, $x, 1) = $newbase;
			$counter++;
		}
	}
	$counter=0;
}

###############################

sub addIndels {

	# target is number of bases to mutate
	while ($counter <= $target ) {

		# pick a random position along sequence 
		# search between 0 < (length - 1) like an index
		$x = round(rand($length-1));

		# 0 = deletion, 1 = insertion
		$y = round(rand(1));
			
		if ($y == 0) { # delete one base
			substr($seq, $x, 1) = "";
			@seq = split(//, $seq);
			$length = scalar(@seq);
		}
		else { # insert one base
			$z = round(rand(3));
			$newbase = $ACGT[$z];
			substr($seq, $x, 1) = $newbase;
			@seq = split(//, $seq);
			$length = scalar(@seq);
		}
		$counter++;
	}
	$counter=0;
}

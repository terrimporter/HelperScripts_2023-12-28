#!/usr/bin/perl
#March 27, 2013 by Terri Porter
#Clean up BOLD sequences from Joel to meet the standards I used for the GenBank training sequences
#usage perl cleanup_COI.plx file.fas

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $kingdom = 'Metazoa';
my $phylum;
my $class;
my $order;
my $family;
my $genus;
#my $species; ### don't worry if no species assignment was available, be lenient here ###
my $gbLine;
my $gb;
my $lineage;
my $seq;
my $length;
my $minlength = 500; ### minimum sequence length for training set ###
my $filename;
my $nextline;
my $j;

#declare array
my @in;
my @line;
my @seq;
my @gbLine;

#decalre hash
my %lineage; #indexed by gb
my %seq; #indexed by gb

open (IN, "<", $ARGV[0]) || die "Error cannot open file.fas: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g;
		@line = split(/\|/, $line);
		$phylum = $line[0];
		if ($phylum eq '_') { #handle missing rank assignments
			$phylum = 'undef'.'_'.$kingdom;
		}
#		print "phylum:$phylum\n";
		$class = $line[1];

		if ($class eq '_') {
			$class = 'undef'.'_'.$phylum;
		}
#		print "class:$class\n";
		$order = $line[2];

		if ($order eq '_') {
			$order = 'undef'.'_'.$class;
		}
#		print "order:$order\n";
		$family = $line[3];

		if ($family eq '_') {
			$family = 'undef'.'_'.$order;
		}
#		print "family:$family\n";
		$genus = $line[4];

		if ($genus eq '_') { #throw this sequence away if no genus level classification provided
			$i++;
			$line=();
			@line=();
			$phylum=();
			$class=();
			$order=();
			$family=();
			$genus=();
			next;
		}
#		print "genus:$genus\n";
		$gbLine = $line[7];
#		print "gbline:$gbLine\n";
		@gbLine = split(/\_/, $gbLine);
		$gb = $gbLine[0];
#		print "gb:$gb\n";

		$lineage = $kingdom.';'.$phylum.';'.$class.';'.$order.';'.$family.';'.$genus;
#		print "lineage:$lineage\n";
		#hash lineage
		$lineage{$gb} = $lineage;

		$j = $i+1;
		$nextline = $in[$j];
		chomp $nextline;

		$seq = $nextline;
#		print "seq:$seq\n";

		if ($seq =~ /N|B|[D-F]|[H-S]|[U-Z]/i) { #missing data or non-nucleotide characters
			delete $lineage{$gb}; #not suitable for the training set so remove
		}

		if ($seq =~ /-/) { #remove all gaps (starting, internal, trailing)
			$seq =~ s/-//g;			
		}
		@seq = split(//, $seq);
		$length = scalar(@seq);
#		print "length:$length\n";

		if ($length > $minlength) {
			$seq{$gb} = $seq;
		}
	}
	elsif ($line =~ /^\s+/) {
		$i++;
		next;
	}
	
	$i++;
	$line=();
	@line=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$gbLine=();
	@gbLine=();
#	$gb=();
	$j=();
}
$i=0;

$filename = $ARGV[0].".parsed";

#report number of sequences with sufficient lineage annotations
while ( ($gb, $lineage) = each(%lineage) ) {
	$i++;
}
print "\nThere were $i accessions with suffucient lineage annotations.\n";
$i=0;

while ( ($gb, $seq) = each(%seq) ) {
	$i++;
}
print "There were $i accessions with sequences that passed quality filters.\n";
$i=0;

open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";

while ( ($gb,$lineage) = each(%lineage) ) {
	
	if (exists $seq{$gb}) {
		$seq = $seq{$gb};
		print OUT ">$gb $lineage\n$seq\n";
		$i++;
	}
	else {
		print "Sequence $gb did not pass quality filters\n";
	}
	$seq=();
}
close OUT;
print "There were $i records in total that passed all filters\n";

#!/usr/bin/perl
#May 1, 2012 by Terri Porter
#Script to exlude kmer.fasta 'baits' with <= 40% GC content.  Start with simple filtering across whole 100bp region.  Consider inreasing sensitivity by using a sliding window approach.
#USAGE perl GC_filter.plx kmer.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $id;
my $j;
my $seq;
my $base;
my $count=0;

#declare array
my @baits;
my @seq;

#declare hash
my %bait;

open (IN, "<", $ARGV[0]) || die "Error cannot read infile: $!\n";
@baits = <IN>;
close IN;

#print "arrayBaits:@baits\n";#test

#add bait seqs to hash
while ($baits[$i]) {
	$line = $baits[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w+)/;
		$id = $1;
		
		$j=$i+1;
		$seq = $baits[$j];
		chomp $seq;

		$bait{$seq} = $id; #this automatically overwrites the value of duplicate keys
	}
	$i++;
	$line=();
	$id=();
	$j=();
	$seq=();
}
$i=0;
$j=0;

open (OUT, ">>", "GC_baits.fasta") || die "Error cannot write to outfile: $!\n";

#filter out baits with <40% GC content
while (my ($seq,$id) = each (%bait)) {
	@seq = split('',$seq);

	while ($seq[$j]) {
		$base = $seq[$j];
		 
		if ($base =~ /(G|C)/i) {
			$count++;
		}
		else {
			$j++;
			next;
		}
		$j++;
		$base=();
	}
	$j=0;

	if ($count >= 40) {#####edit GC content cutoff here, assumes 100 bp baits#####
		print OUT ">$id\n$seq\n";
	}
	$count=0;
	@seq=();
}

close OUT;

#!/usr/bin/perl
#May 5, 2013 edited to also remove A/T runs of 8bp +
#May 1, 2012 by Terri Porter
#Script to exlude kmer.fasta 'baits' with <= 40% GC content.  Start with simple filtering across whole 100bp region.  Consider inreasing sensitivity by using a sliding window approach.  Automatically removes duplciates.
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
my $GCcutoff = 40; ### reset minimum GC content filter here, assumes 100bp baits ###

#declare array
my @baits;
my @seq;

#declare hash
my %bait; #indexed by sequence to automatically remove/overwrite duplicates
my %gc40; #indexed by id
#my %gc40noruns; #indexed by id

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

open (OUT, ">>", "GC_40plus.fasta") || die "Error cannot open first outfile: $!\n";

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

	if ($count >= $GCcutoff ) {
		print OUT ">$id\n$seq\n";
		$gc40{$id} = $seq;
	}
	$count=0;
	@seq=();
}

close OUT;

#May 7, 2013 add another filtering step
#Remove baits with GC stretches 8bp+
open (OUT2, ">>", "GC_40plus_noATruns.fasta") || die "Error cannot open second outfile: $!\n";

while ( ($id, $seq) = each (%gc40) ) {
	if ($seq =~ /(A|T){8,}/) { #look for runs of 8 or more A/T bases
		print "found a run: $seq\n";
		next;
	}
	else {
		print OUT2 ">$id\n$seq\n";
#		$gb40noruns{$id} = $seq;
	}
}
close OUT2;

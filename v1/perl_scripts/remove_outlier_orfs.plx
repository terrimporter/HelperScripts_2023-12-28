#!/home/terri/miniconda3/envs/myenv/bin/perl
# Dec. 12, 2019 by Teresita M. Porter
# Tweak to account for fasta headers formatted like: Amplicon_OTUid
# Script to grab longest CDS from ORFfinder nucleotide output (option 1)
# only keep CDS with lengths within 25/75th percentile -/+ 1.5*IQR
# print to STDOUT here, in snakemake, redirect to an outfile
# already working with longest OTUs, don't bother sorting by size
# USAGE perl get_longest_orf.plx cds.fasta.tmp limits.txt

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $ampliconOtu;
my $j;
my $seq;
my $length;

my $lower;
my $upper;

# declare array
my @in;
my @seq;
my @range;

# declare hash
my %seq; #key1 = ampliconOtu, key2 = length, value = seq

open (IN, "<", $ARGV[0]) || die "Can't open infile: $!\n";
@in = <IN>;
close IN;

# parse file to get fasta header stats first
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~/^>/) {
		$line =~ s/^>//g; # >Amplicon_OTUid
		$ampliconOtu = $line;
		$j = $i+1;
		$seq = $in[$j];
		chomp $seq;
		@seq = split(//,$seq);
		$length = scalar(@seq);
		$seq{$ampliconOtu}{$length} = $seq; #create hash of hashes
		$i++;
		next;
	}
	else {
		$i++;
		next;
	}
}
$i=0;
$j=();
		

# loop through each otu, get orf id for longest length, print otu and seq to outfile

# get 25/75th percentile -/+ 1.5*IQR range
open (IN2, "<", $ARGV[1]) || die "Can't open limits.txt:$!\n";
@range = <IN2>;
close IN2;

$lower = $range[0];
chomp $lower;
$upper = $range[1];
chomp $upper;

# print ORFs within range (exclude outlier orfs)
foreach $ampliconOtu (keys %seq) {
	foreach $length (keys %{$seq{$ampliconOtu}}) {
			
		if ($length >= $lower && $length <= $upper) {
			$seq = $seq{$ampliconOtu}{$length};
			print STDOUT ">$ampliconOtu\n$seq\n";
			$i++;
			next;
		}
	}
}

#!/usr/bin/perl
# Teresita Porter, Sept. 13, 2020
# Script to add QIIME formatted lineage back to the headers that disappeared after 99% clustering
# USAGE perl map_lineage_to_clusters.plx GenBank_BOLD.fasta GenBank_BOLD.99.fasta

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $accession;
my $lineage;
my $outfile = "GenBank_BOLD.99.mapped.fasta";
my $flag=0;
my $seq;
my $catseq;

# declare array
my @map;
my @line;
my @clusters;

# declar hash
my %map; #key = accession, value = lineage

open (IN, "<", $ARGV[0]) || die "Error cannot open first infile: $!\n";
@map = <IN>;
close IN;

# hash the data from mapping FASTA file
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g;
		@line = split(' ',$line);
		$accession = $line[0];
		$lineage = $line[1];
		$map{$accession}=$lineage;
	}
	$i++;

}
$i=0;

# create new outfile with lineages added
open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

open (IN2, "<", $ARGV[1]) || die "Error cannot open second infile: $!\n";
@clusters = <IN2>;
close IN2;

# be sure to convert to a STRICT FASTA format with sequence on ONE line
while ($clusters[$i]) {
	$line = $clusters[$i];
	chomp $line;

	if ($line =~ /^>/) {
		# print the previous full sequence on one line
		if ($flag==1) {
			print OUT $seq."\n";
			$seq="";
			
		}
		# map the lineage and print it
		$line =~ s/^>//g;
	
		if (exists $map{$line}) {
			$lineage = $map{$line};
			print OUT ">".$line." ".$lineage."\n";
		}
		else {
			print "Can't find accession $line\n";
		}

		$flag = 1;

	}
	elsif ($flag=1) {
		$catseq = $seq.$line;
		$seq = $catseq;
	}

	$i++;
}
$i=0;

#print the final seq
print OUT $seq."\n";
close OUT;

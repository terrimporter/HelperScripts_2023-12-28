#!/usr/bin/perl
# Teresita M. Porter, May 3/21
# Script to add lineage back to FASTA from vsearch, also get rid of line breaks in sequence
# skip over accessions where lineage contains sp. or cf. or aff.
# also remove starting and trailing N's in the BOLD seqs (dashes were removed earlier on)
# screen out any BOLD seqs < 500 bp
# also create taxid.parsed
# USAGE perl add_lineage_to_fasta.plx bat.fasta.derep mapping

use strict;
use warnings;
use Data::Dumper;

# vars
my $i = 0;
my $line;
my $header;
my $acc;
my $lineage;
my $seq;
my $outfile = "testNBC.fasta";
my $outfile2 = "taxid.parsed";
my $flag=0;
my $len;

# arrays
my @fas;
my @map;
my @line;
my @header;

# hashes
my %map; # key = acc, value = lineage
my %lineage; #key = lineage, value = 1

open (IN, $ARGV[0]) || die "Can't open infile1: $!\n";
@fas = <IN>;

open (MAP, $ARGV[1]) || die "Can't open mapping file: $!\n";
@map = <MAP>;

# hash the mappinng file
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$acc = $line[0];
	$lineage = $line[1];

	if ($lineage =~ /cf\./) {
		$i++;
		next;
	}
	elsif ($lineage =~ /sp\./) {
		$i++;
		next;
	}
	elsif ($lineage =~ /aff\./){
		$i++;
		next;
	}
	elsif (! defined $lineage) { # check if lineage is defined
		$i++;
		next;
	}
	else {
		$map{$acc} = $lineage;
	}
	$i++;
	$acc=();
	$lineage=();
}
$i=0;

# testing
#print Dumper(\%map);

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";
open (OUT2, ">>", $outfile2) || die "Error cannot open outfile2: $!\n";

while ($fas[$i]) {
	$line = $fas[$i];
	chomp $line;

	if ($line =~ /^>/ && $flag==0) { # the first header or the first header after skipping a header
		$acc = $line;
		$acc =~ s/^>//g;
		if (exists $map{$acc}) {
			$lineage = $map{$acc};		
			$header = ">$acc\t$lineage";
			$seq="";
			$flag=1;
		}
		else {
			print "Can't find lineage for $acc\n";
			$seq="";
			$flag=0;
		}
	}
	elsif ($line =~ /^>/ && $flag==1) { # the rest of the headers
		chomp $header;
		# remove starting and trailing N's from BOLD seqs
		$seq =~ s/^-*//g;
		$seq =~ s/-*$//g;
		$len = length $seq;

		if ($len >= 500) { # avoid short seqs
			print OUT $header."\n".$seq."\n"; # print previous entry to fasta file
			@header = split(/\t/, $header);
			$lineage = $header[1];
			$lineage =~ s/;/\t/g;

			if (!exists $lineage{$lineage}) {
				$lineage{$lineage} = 1;
				print OUT2 $lineage."\n"; # print lineage to taxid.parsed
			}
		}

		$acc = $line;
		$acc =~ s/^>//g;
		if (exists $map{$acc}){
			$lineage = $map{$acc};
			$header = ">$acc\t$lineage";
			$seq="";
			$flag=1;
		}
		else {
			print "Can't find lineage for $acc\n";
			$seq="";
			$flag=0;
		}
	}
	elsif ($flag==1) { # only process seq if good acc
		$seq = $seq.$line;
	}
	$i++;	

}
$i=0;

# print last record
if ($flag==1) {
	print OUT $header."\n".$seq."\n";
	@header = split(/\t/, $header);
	$lineage = $header[1];
	$lineage =~ s/;/\t/g;

	if (!exists $lineage{$lineage}) {
		print OUT2 $lineage."\n"; # print lineage to taxid.parsed
	}
}

close OUT;
close OUT2;

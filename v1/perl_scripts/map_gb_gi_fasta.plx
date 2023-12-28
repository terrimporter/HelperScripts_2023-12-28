#!/usr/bin/perl
#March 16,2011 by Terri Porter
#Script to take a map file of gb accessions and corresponding gi numbers, and add a gi number to the matching gb accession in the fasta file
#usage map_gb_gi_fasta.plx gb_gi.query ITS.fasta


use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb_acc;
my $j=0;
my $line2;
my $gi;

#declare array
my @map;
my @fasta;
my @line;
my @line2;


open (MAP, "<", $ARGV[0]) || die ("Error cannot read map file: $!\n");
@map = <MAP>;
close MAP;

open (FASTA,"<",$ARGV[1]) || die ("Error cannot read fasta file: $!\n");
@fasta = <FASTA>;
close FASTA;

open (OUT,">>","ITS.fasta.mapped") || die ("Error cannot write mapped fasta file: $!\n");

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ s/>//;
		@line = split(/\|/,$line);
		$gb_acc = $line[0];

		while ($map[$j]) {
			$line2 = $map[$j];
			chomp $line2;
			if ($line2 =~ /$gb_acc/) {
				@line2 = split(/\t/,$line2);
				$gi = $line2[1];
				$gi =~ s/GI://;
				print OUT ">$gi|$gb_acc\n";
			}
			$j++;
		}
		$j=0;
	}
	else {
		print OUT "$line\n";
	}
	$i++;
}
close OUT;

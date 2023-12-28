#!/usr/bin/perl
#Feb. 15, 2013 by Terri Porter
#Script to grab fasta headers from *.muscle and put these back on *.trimal (trimal stripped everything after the first space!
#usage perl map_muscle_to_trimal_header.plx file.muscle file.trimal

use strict;
use warnings;

#declare var
my $i=0;
my $muscleline;
my $numseq=0;
my $trimalline;
my $numseq2=0;

#declare array
my @muscle;
my @trimal;

#declare hash
my %muscle;

open (MUSCLE, "<", $ARGV[0]) || die "Error cannot open muscle file: $!\n";
@muscle = <MUSCLE>;
close MUSCLE;

open (TRIMAL, "<", $ARGV[1]) || die "Error cannot open trimal file: $!\n";
@trimal = <TRIMAL>;
close TRIMAL;

#hash muscle headers
while ($muscle[$i]) {
	$muscleline = $muscle[$i];
	chomp $muscleline;

	if ($muscleline =~ /^>/) {
		$numseq++;
		$muscle{$numseq} = $muscleline;
	}
	$i++;
}
$i=0;

open (OUT, ">>", "trimal.fixed") || die "Error cannot open outfile: $!\n";

#reprint trimal with headers from muscle
while ($trimal[$i]) {
	$trimalline = $trimal[$i];
	chomp $trimalline;

	if ($trimalline =~ /^>/) {
		$numseq2++;
		$muscleline = $muscle{$numseq2};
		print OUT "$muscleline\n";
	}
	else {
		print OUT "$trimalline\n";
	}
	$i++;
}
$i=0;

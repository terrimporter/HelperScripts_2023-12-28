#!/usr/bin/perl
#Nov. 25, 2011 by Terri Porter
#Script to get features for 33 parent sequences indexed by GI from features file indexed by GB
#usage perl get_features.plx features.txt.parsed gb.query gi.query overlap.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $gb_v;
my $gb;
my $gi;
my $gensp;
my $line;
my $lineage;
my $header;

#declare array
my @feat;
my @gb_v;
my @gi;
my @fasta;
my @gb;
my @line;

#declare hash
my %gi_gb;
my %gb_lineage;
my %gb_gensp;

open (FEAT,"<",$ARGV[0]) || die "Error cannot open features.txt.parsed: $!\n";
@feat=<FEAT>;
close FEAT;

open (GB,"<",$ARGV[1]) || die "Error cannot open gb.query: $!\n";
@gb_v=<GB>;
close GB;

open (GI,"<",$ARGV[2]) || die "Error cannot open gi.query: $!\n";
@gi=<GI>;
close GI;

open (FASTA,"<",$ARGV[3]) || die "Error cannot open overlap.fasta: $!\n";
@fasta=<FASTA>;
close FASTA;

#remove versions from gb
while ($gb_v[$i]) {
	$gb_v = $gb_v[$i];
	chomp $gb_v;
	$gb_v =~ s/\.\d+$//;
	push(@gb,$gb_v);
	$i++;
}
$i=0;

#create gi_gb hash map
while ($gb[$i]) {
	$gb=$gb[$i];
	$gi=$gi[$i];
	chomp $gi;
	$gi_gb{$gi} = $gb;
	$i++;
}
$i=0;

#create gb_lineage and gb_gensp hash maps
while ($feat[$i]) {
	$line = $feat[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gb = $line[0];
	$gensp = $line[1];
	$lineage = $line[5];
	$gb_gensp{$gb} = $gensp;
	$gb_lineage{$gb} = $lineage;
	$i++;
	@line=();
	$gb=();
	$gensp=();
	$lineage=();
}
$i=0;

open (OUT,">>","parsed_features.txt") || die "Error cannot write to parsed_features.txt: $!\n";

#get lineage and gensp for each gi in fasta file
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/^>//;
		$gi = $header;
		$gb = $gi_gb{$gi};
		$lineage = $gb_lineage{$gb};
		$gensp = $gb_gensp{$gb};
		print OUT "$gi\t$gb\t$lineage\t$gensp\n";
		$gi=();
		$gb=();
		$lineage=();
		$gensp=();
	}
	$i++;
}
$i=0;
close OUT;	

#!/usr/bin/perl
#Nov. 27, 2013 by Terri Porter
#Script to grab gi's for each gb in a fasta file header
#USAGE perl grab_gi_for_each_gb_testNBC.plx gb_gi.map.uniq testNBC.fasta

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $gb;
my $gi;

#declare array
my @map;
my @fasta;
my @line;

#declare hash
my %map; #key = gb, value = gi

open (MAP, "<", $ARGV[0]) || die "Error cannot open map file: $!\n";
@map = <MAP>;
close MAP;

open (FASTA, "<", $ARGV[1]) || die "Error cannot open fasta file: $!\n";
@fasta = <FASTA>;
close FASTA;

#hash map file
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$gb = $line[0];
	$gi = $line[1];
	$map{$gb} = $gi;
	
	$i++;
	$line=();
	@line=();
	$gb=();
	$gi=();
}
$i=0;

#parse fasta

open (OUT, ">>", "testNBC.fasta.gi") || die "Error cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w+)\s{1}/;
		$gb = $1;
		
		if (exists $map{$gb}) {
			$gi = $map{$gb};
			print OUT $gi."\n";
		}
		else {
			print "Cannot find gb $gb\n";
		}
	}

	$i++;
	$line=();
	$gb=();
	$gb=();

}
$i=0;
close OUT;

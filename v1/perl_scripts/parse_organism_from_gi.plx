#!/usr/bin/perl
#April 7, 2011 by Terri Porter
#Script to parse hit_gb.txt from get_gi_and_lineage_for_each_gb.plx and grab just gb and organism name
#usage $perl parse_organism_from_gi.plx hit_gb.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $hit_gb;
my $hit_genus;
my $hit_species;
my $flag=1;

#declare array
my @in;

open (IN,"<",$ARGV[0]) || die ("Error cannot read from hit_gb.txti: $!\n");
@in = <IN>;
close IN;

open (OUT,">>","hit_gb.txt.organism") || die ("Error cannot write to hit_gb.txt.organism: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	if ($line =~ /LOCUS/) {
		$flag=1;
	}
	elsif ($flag==1) {
	if ($line =~ /VERSION/) {
		$line =~ /VERSION\s+(\w+)\.\d+/;
		$hit_gb = $1;
	}
	elsif ($line =~ /ORGANISM/) {
		$line =~ /ORGANISM\s+(\w+)\s+(\w+)/;
		$hit_genus = $1;
		$hit_species = $2;
		print OUT "$hit_gb\t$hit_genus\t$hit_species\n";
		$hit_gb=();
		$hit_genus=();
		$hit_species=();
		$flag=0;
	}
	}
	$i++;
}
close OUT;

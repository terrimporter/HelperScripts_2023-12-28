#!/usr/bin/perl
#Sept. 22, 2011 by Terri Porter
#also removes lines that begin with "nil" (these were duplicate strains/isolates, only the longest seq was kept)
#Script to parse features.txt with too many entires using gb.list originally created by LSU_entrez.plx
#also remove accessions starting with NG_ because these are refseqs annotated by GenBank staff that are duplicates of sequences I already have
#usage perl parse_features_with_gb_list.plx gb.query features.txt

use strict;
use warnings;

#declare var
my $i=0;
my $j=0;
my $gb;
my $gb_short;
my $line;

#declare array
my @gb;
my @feat;
my @gb_split;

#declare hash

open (GB,"<",$ARGV[0]) || die "Error cannot open gb.query: $!\n";
@gb = <GB>;
close GB;
print "@gb\n";#test

open (FEAT,"<",$ARGV[1]) || die "Error cannot open features.txt: $!\n";
@feat = <FEAT>;
close FEAT;

open (OUT,">>","features.txt.parsed") || die "Error cannot write to features.txt.parsed: $!\n";

while ($gb[$i]) {
	$gb = $gb[$i];
	chomp $gb;
	$gb =~ /(\w+)\.\d+/;
	$gb_short = $1;
	print "gb_short:$gb_short\n";#test
	#	print "array:@gb_split\n";#test
#	$gb_short = $gb_split[0];
#	print "gb_short:$gb_short\n";#test

	while ($feat[$j]) {
		$line = $feat[$j];
		chomp $line;
		if ($line !~ /^NG_/) {
			if ($line =~ /^$gb_short/) {
				print OUT "$line\n";
			}
		}
		$j++;
	}
	$j=0;
	$i++;
}
close OUT;

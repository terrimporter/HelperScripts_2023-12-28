#!/usr/bin/perl
#May 16, 2011 by Terri Porter
#Script to parse out gi number from megan_readid_leaf.txt file
#usage $perl parse_taxonid_from_megan_taxonid.plx megan_readid_taxonid.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $taxonid;

#declare array
my @megan;
my @line;

open (IN,"<",$ARGV[0]) || die ("Error cannot read megan_readid_taxonid.txt: $!\n");
@megan = <IN>;
close IN;

open (OUT,">>","megan.taxonid") || die ("Error cannot write to megan.gi: $!\n");

while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;

	if ($line =~ /^\d+/) {
		@line = split(/,/,$line);
		$taxonid = $line[1];
		$taxonid =~ s/ //;
		print OUT "$taxonid\n";
	}
	$i++;
}

close OUT;

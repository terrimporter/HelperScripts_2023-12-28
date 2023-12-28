#!/usr/bin/perl
#June 13, 2012 by Terri Porter
#Script to fix NC_# and change to NC# so that compare_taxonomy_maps.plx counts everything right
#usage perl fix_NC_gb_org.plx gb_org.map

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $new;

#declare array
my @in;

open (IN, "<", $ARGV[0]) || die "Error cannot open gb_org.map infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "gb_org_NCfixed.map") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^NC/) {
		$line =~ /^(NC)_(\w+)/;
		$new = $1.$2;
		print OUT "$new\n";
	}
	else {
		print OUT $line."\n";
	}

	$i++;
	$line=();
	$new=();
}
$i=0;
close OUT;

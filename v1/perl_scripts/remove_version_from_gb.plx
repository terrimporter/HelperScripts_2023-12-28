#!/usr/bin/perl
#April 7, 2011 by Terri Porter
#Script to remove version part of gb from a list of gb
#usage $perl remove_version_from_gb.plx hit_gb.list

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $new_gb;

#declare array
my @gb;
my @line;

open (GB,"<",$ARGV[0]) || die ("Error cannot read from hit_gb.list: $!\n");
@gb=<GB>;
close GB;

open (OUT,">>","hit_gb.list.noversion") || die ("Error cannot write to hit_gb.list.noversion: $!\n");

while ($gb[$i]) {
	$line = $gb[$i];
	chomp $line;

	@line = split(/\./,$line);
	$new_gb = $line[0];

	print OUT "$new_gb\n";
	$i++;
}
close OUT;

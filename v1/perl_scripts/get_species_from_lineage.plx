#!/usr/bin/perl
#June 7, 2011 by Terri Porter
#Script to grab species name only from query_lineage.txt
#usage perl get_species_from_lineage.plx hit_lineage.txt.mapped

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $species;
my $map;

#declare array
my @query;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@query = <IN>;
close IN;

while ($query[$i]) {
	$line = $query[$i];
	chomp $line;
	if ($line =~ /^\d+/) {
		$line =~ /^(\d+)/;
		$map = $1;
		if ($line =~ /\w+\.$/) {
			$line =~ /(\w+)\.$/;
			$species = $1;
			print "$map\t$species\n";
		}
		elsif ($line =~ /\w+\.\.$/) {
			$line =~ /(\w+)\.\.$/;
			$species = $1;
			print "$map\t$species\n";
		}
	}
	$i++;
}
close IN;

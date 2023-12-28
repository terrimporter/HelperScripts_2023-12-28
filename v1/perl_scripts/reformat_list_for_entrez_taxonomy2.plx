#!/usr/bin/perl
#May 9, 2013 edited to reformat Family into a form that can be used with ebot_taxonomy.plx
#May 14, 2012 by Terri Porter
#Script to reformat Genus\tSpecies into a form that can be used with ebot_taxonomy.plx
#usage perl reformat_list_for_ebot_taxonomy.plx page.species

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $genus;
my $species;
my $family;
my $max; #total number of species
my $target; #index number of last species in array
my $filename; #generate custom filename

#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0].".txt";

open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";

$max = scalar(@in);
$target = $max-1;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$family = $line[0];
#	$species = $line[1];

	if ($i < $target) {
		print OUT "\"$family\"[Organism] OR ";
	}
	elsif ($i == $target) {
		print OUT "\"$family\"[Organism]";
	}

	$i++;
	$line=();
	@line=();
	$family=();
#	$species=();
}
$i=0;

close OUT;

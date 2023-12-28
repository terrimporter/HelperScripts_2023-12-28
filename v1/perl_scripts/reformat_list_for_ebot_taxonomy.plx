#!/usr/bin/perl
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

#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "list_for_ebot.txt") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$genus = $line[0];
	$species = $line[1];
	print OUT "%22$genus $species%22[Organism]+OR+";

	$i++;
	$line=();
	@line=();
	$genus=();
	$species=();
}
$i=0;

close OUT;

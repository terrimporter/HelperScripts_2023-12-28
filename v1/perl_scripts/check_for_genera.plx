#!/usr/bin/perl
# Jan. 10, 2019 by Terri M. Porter
# Script to search for taxa in the COI Classifier reference database.  Check the taxonomy file because its a lot smaller than the FASTA file.
# USAGE perl check_for_genera.plx targetGenera.txt testNBC.taxonomy

use strict;
use warnings;

# declare var
my $taxon;
my $line;
my $i=0;
my $species;
my $outfile = "GenusSearch.results";

# declare array
my @target;
my @ref;
my @line;

# declare hash
my %target; # key = taxon, value = 1
my %ref; # key = taxon, value = 1

# open infiles
open (IN, "<", $ARGV[0]) || die "Cannot open first infile: $!\n";
@target = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Cannot open second infile: $!\n";
@ref = <IN2>;
close IN2;

# hash the target taxa (Genus)
foreach $taxon (@target) {
	chomp $taxon;
	$target{$taxon} = 1;
}

# hash reference taxonomy
while ($ref[$i]) {
	$line = $ref[$i];
	chomp $line;

	if ($line =~ /genus$/) {
		@line = split(/\*/,$line);
		$species = $line[1];
		$ref{$species} = 1;
		$i++;
		next;
	}
	else {
		$i++;
		next;
	}
}
$i=0;

# check for target taxa in reference taxonomy

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

while (my ($key, $value) = each(%target)) {
	if (exists $ref{$key}) {
		print OUT "$key\tPRESENT\n";
	}
	else {
		print OUT "$key\tMISSING\n";
		next;
	}

}
close OUT;

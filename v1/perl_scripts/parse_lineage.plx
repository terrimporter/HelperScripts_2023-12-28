#!/usr/bin/perl
#March 22, 2011 by Terri Porter
#Script to parse GenBank lineage info for each gi, find family, order, and class level info
#usage $perl parse_lineage.plx parsed_genbankrecords.outfile

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $j=0;
my $rank;
my $family;
my $suborder;
my $order;
my $subclass;
my $class;
my $genus_gi;
my $gi;
my $genus;
my $subphylum;
my $phylum;
my $subkingdom;
my $kingdom;
my $superkingdom;

#declare arrays
my @lineage;
my @line;
my @genus_gi;

open (IN,"<",$ARGV[0]) || die ("Error cannot read lineage info: $!\n");
@lineage = <IN>;
close IN;

open (OUT, ">>", "parsed_lineage.txt") || die ("Error cannot write parsed lineage info: $!\n");

while ($lineage[$i]) {
	$line = $lineage[$i];
	chomp $line;
	@line = split(/;/,$line);
	
	while ($line[$j]) {
		$rank = $line[$j];

		if ($rank =~ /aceae$/) {
			$family = $rank;
		}
		elsif ($rank =~ /ineae$/) {
			$suborder = $rank;
		}
		elsif ($rank =~ /ales$/) {
			$order = $rank;
		}
		elsif ($rank =~ /mycetidae$/) {
			$subclass = $rank;
		}
		elsif ($rank =~ /mycetes$/) {
			$class = $rank;
		}
		elsif ($rank =~ /\d+$/) {
			$genus_gi = $rank;
			@genus_gi = split(/\t/,$genus_gi);
			$gi = $genus_gi[1];
			$genus = $genus_gi[0];
		}
		elsif ($rank =~ /mycotina$/) {
			$subphylum = $rank;
		}
		elsif ($rank =~ /mycota$/) {
			$phylum = $rank;
		}
		elsif ($rank =~ /Dikarya$/) {
			$subkingdom = $rank;
		}
		elsif ($rank =~ /Fungi$/) {
			$kingdom = $rank;
		}
		elsif ($rank =~ /(Eukaryota$)/) {
			$superkingdom = $rank;
		}
		$j++;
	}
	$j=0;

	print OUT "$gi\t$genus\t$family\t$suborder\t$order\t$subclass\t$class\t$subphylum\t$phylum\t$subkingdom\t$kingdom\t$superkingdom\n";
	$i++;
	$gi=undef;
	$genus=undef;
	$family=undef;
	$suborder=undef;
	$order=undef;
	$subclass=undef;
	$class=undef;
	$subphylum=undef;
	$phylum=undef;
	$subkingdom=undef;
	$kingdom=undef;
	$superkingdom=undef;
												
}

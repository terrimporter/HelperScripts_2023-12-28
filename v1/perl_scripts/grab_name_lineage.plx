#!/usr/bin/perl
#Sept.26,2011 by Terri Porter
#for each accession in primer_400.fasta grab name and lineage from features.txt.parsed
#usage perl grab_name_lineage.plx primer_400.fasta features.txt.parsed

use strict;
use warnings;

#declare var
my $i=0;
my $id;
my $j=0;
my $k=0;
my $line;
my $organism;
my $lineage;

#declare array
my @id;
my @features;
my @ids_to_search;
my @line;

open (ID,"<",$ARGV[0]) || die "Error cannot open primer_400.fasta: $!\n";
@id = <ID>;
close ID;

open (FEATURES,"<",$ARGV[1]) || die "Error cannot open features.txt.parsed: $!\n";
@features = <FEATURES>;
close FEATURES;

open (OUT,">>","id_lineage.txt") || die "Error cannot write to id_lineage.txt: $!\n";

while($id[$i]) {
	$id = $id[$i];
	chomp $id;

	if ($id =~ /^>/) {
		$id =~/^>(\w+)/;
		$id = $1;
		push(@ids_to_search,$id);
	}
	$i++;
}

while ($ids_to_search[$j]) {
	$id = $ids_to_search[$j];

	while ($features[$k]) {
		$line = $features[$k];
	
		if ($line =~ /^$id/) {
			@line = split(/\t/,$line);
			$organism = $line[1];
			$lineage = $line[5];
			print OUT "$id\t$organism\t$lineage\n";
		}
		$k++;
	}
	$k=0;
	$j++;
}
close OUT;			

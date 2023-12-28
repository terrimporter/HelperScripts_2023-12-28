#!/usr/bin/perl
#Nov. 23, 2016 by Terri Porter, script to check a list of taxa against the RDP classifier training set for inclusion
#USAGE perl xcheck_rdp.plx genera.txt testNBC.fasta

use strict;
use warnings;

#declar var
my $i=0;
my $j=0;
my $line;
my $field1;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $gb;
my $superkingdom;

#declare array
my @genera;
my @nbc;
my @line;
my @field1;

#declare hash
my %nbc; #keys are genera for easy searching

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@genera = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@nbc = <IN2>;
close IN2;

while ($nbc[$j]) {
	$line = $nbc[$j];
	chomp $line;

	if ($line =~ /^>/) {

		@line = split(";",$line);
		$field1 = $line[0];
		$field1 =~ s/>//g;
#		print $field1."\n";

		$kingdom = $line[1];
		$phylum = $line[2];
		$class = $line[3];
		$order = $line[4];
		$family = $line[5];
		$genus = $line[6];
#		print $genus."\n";#test

		@field1 = split(" ",$field1);
		$gb = $field1[0];
		$superkingdom = $field1[1];

		$nbc{$genus} = "1";

	}

	$j++;

}
$j=0;

open (OUT1, ">>", "genera_in_nbc.txt") || die "Cannot open outfile1: $!\n";
open (OUT2, ">>", "genera_not_in_nbc.txt") || die "Cannot open outfile2: $!\n";

while ($genera[$i]) {
	$line = $genera[$i];
	chomp $line;

	if (exists $nbc{$line}) {
		print OUT1 "$line\n";
	}
	else {
		print OUT2 "$line\n";
	}
	$i++;
}
$i=0;
close OUT1;
close OUT2;

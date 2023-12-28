#!/usr/bin/perl
#Oct. 3, 2013 Script to parse NBC classifier output, specifically for seqs classified using Insecta COI trained classifier
#USAGE perl parse_NBC.plx file.out

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $id;
my $kingdom;
my $kingdom_bp;
my $phylum;
my $phylum_bp;
my $class;
my $class_bp;
my $order;
my $order_bp;
my $family;
my $family_bp;
my $genus;
my $genus_bp;
my $cutoff=0.7; ### adjust genus rank bootstrap cutoff here ###


#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open NBC file: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "NBC.parsed") || die "Error cannot open NBC parsed file: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if (length($line) > 0) {
		@line = split(/\t/,$line);
		$id = $line[0];
		$kingdom = $line[2];
		$kingdom_bp = $line[4];
		$phylum = $line[5];
		$phylum_bp = $line[7];
		$class = $line[8];
		$class_bp = $line[10];
		$order = $line[11];
		$order_bp = $line[13];
		$family = $line[14];
		$family_bp = $line[16];
		$genus = $line[17];
		$genus_bp = $line[19];

		if ($genus_bp >= $cutoff) {
			print OUT "$id\t$kingdom\t$kingdom_bp\t$phylum\t$phylum_bp\t$class\t$class_bp\t$order\t$order_bp\t$family\t$family_bp\t$genus\t$genus_bp\n";
		}
	}
	$i++;
	$line=();
	@line=();
	$kingdom=();
	$kingdom_bp=();
	$phylum=();
	$phylum_bp=();
	$class=();
	$class_bp=();
	$order=();
	$order_bp=();
	$family=();
	$family_bp=();
	$genus=();
	$genus_bp=();
}
$i=0;
close OUT;

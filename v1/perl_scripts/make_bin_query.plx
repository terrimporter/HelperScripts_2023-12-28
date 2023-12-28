#!/usr/bin/perl
#May 9, 2012 by Terri Porter
#Script to parse Insecta.txt (from parse_tsv.plx) and construct a BIN query at http://www.boldsystems.org/index.php/Public_BarcodeIndexNumber_Home
#usage perl make_bin_query.plx Insecta.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $header;
my $entry;
my $species_reg;
my $BIN;

#declare array
my @in;
my @entry;
my @species_reg;

#declare hash
my %bin;

open (IN,"<",$ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT,">>","BIN.query") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~/^processid/) {
		$header = $line;
		$i++;
		next;
	}
	else {
		$entry = $line;
		@entry = split(/\t/,$entry);
		$species_reg = $entry[6];
		@species_reg = split(" ",$species_reg);
		$BIN = pop(@species_reg);
		$bin{$BIN} = 1;
	}
	$i++;
	$line=();
	$header=();
	$entry=();
	@entry=();
	$species_reg=();
	@species_reg=();
	$BIN=();
}
$i=0;

while(my($BIN,$value) = each(%bin)) {
	print OUT "$BIN ";
}
close OUT;

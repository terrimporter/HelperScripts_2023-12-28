#!/usr/bin/perl
#March 1,2011 by Terri Porter
#first create gb_gi.query $grep "VERSION" myseqs.gb | gawk '{print $2"\t"$3} > gb_gi.query
#Script to parse ITS.fasta, for each gb base name, get corresponding gi from gb_gi.query and make comma-delimited string for blastn -negative_gilist
#usage $perl parse_gb_gi.plx ITS.fasta gb_gi.query

use strict;
use warnings;

#var
my $i=0;
my $line;
my $gb_base;
my $j=0;
my $k=0;
my $base;
my $to_match;
my $gi_part;
my $comma_delim;
my $number;

#array
my @fasta;
my @map;
my @line;
my @gb_base;
my @to_match;
my @gi;

open (FASTA,"<",$ARGV[0]) || die ("Error can't read fasta file:$!|n");
@fasta = <FASTA>;
close FASTA;

open (MAP,"<",$ARGV[1]) || die ("Error can't read gb_gi map file: $!\n");
@map = <MAP>;
close MAP;

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$line =~s/>//;
		@line = split(/\|/,$line);
		$gb_base = $line[0];
		push(@gb_base,$gb_base);
	}
	$i++;
}

#print "array gb_base: @gb_base\n";

while ($gb_base[$j]) {
	$base = $gb_base[$j];

	while ($map[$k]) {
		$to_match = $map[$k];
		chomp $to_match;

		if ($to_match =~ /$base/) {
			@to_match = split(/\t/,$to_match);
			$gi_part = $to_match[1];
			$gi_part =~ s/GI://;
			push(@gi,$gi_part);
		}
		$k++;
	}
	$k=0;
	$j++;
}

#print "array gi: @gi\n";

$comma_delim = join(" ",@gi);
print "$comma_delim\n";
$number = scalar(@gi);
#print "$number elements in gi array\n";

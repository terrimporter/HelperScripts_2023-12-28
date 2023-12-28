#!/usr/bin/perl
# Terri M. Porter, Nov. 23/20
# Use a comma separated mapping file to edit leaf names
# USAGE perl replace_newick_leaf_names.plx map.txt filename.treefile

use strict;
use warnings;
use Data::Dumper;

# scalar
my $i=0;
my $key;
my $value;
my $tempfile = $ARGV[1].".temp";
my $outfile = $ARGV[1].".map";
my $line;
my $newline = "";

# array
my @map;
my @line;
my @tree;
my @in;

# hash
my %map; # key = accession, value = name

open (MAP, "<", $ARGV[0]) || die "Error can't open mapping file:$!\n";
@map = <MAP>;
close MAP;

# hash array
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/,/,$line);
	$key = $line[0];
	if ($key =~ /BOLD:/) {
		$key =~ s/BOLD://;
		print "replacement done\n";
	}
	$line[1] =~ s/ /_/g;
	$value = $key."_".$line[1];
	$map{$key} = $value;
	$i++;
}
$i=0;

# read in the treefile
open (TREE, "<", $ARGV[1]) || die "Error can't open mapping file: $!\n";
@tree = <TREE>;
close TREE;

open (TMP, ">>", $tempfile) || die "Error can't open outfile: $!\n";

# combine multi line newick into a single line
while ($tree[$i]) {
	$line = $tree[$i];
	chomp $line;

	print TMP $line;
	$i++;
}
$i=0;
close TMP;

open (IN, "<", $tempfile) || die "Error can't open tempfile: $!\n";
@in = <IN>;
close IN;

# change names in the treefile
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	foreach $key (keys %map) {
		$value = $map{$key};

		if (length($newline)) {
			$newline =~ s/$key/$value/g;
		}
		else {
			$line =~ s/$key/$value/g;
		}
	}

	$newline = $newline.$line;
	$i++;

}
$i=0;

open (OUT, ">>", $outfile) || die "Error can't open outfile: $!\n";
print OUT $newline."\n";
close OUT;

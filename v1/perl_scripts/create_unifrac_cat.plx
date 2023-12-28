#!/usr/bin/perl
#Nov. 15, 2013 by Terri Porter
#Script to create a category mapping file for fast unifrac
#need env.txt from create_unifrac_map.plx
#usage perl create_unifrac_cat.plx env.txt

use strict;
use warnings;

#declare var
my $line;
my $sample;
my $env;

#declare array
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
open (OUT, ">>", "category.txt") || die "Error cannot open outfile: $!\n";

while (<IN>) {
	$line = $_;
	chomp $line;

	@line = split(/\t/,$line);
	$sample = $line[1];
	$sample =~ s/ //g; #remove any spaces

	if ($sample =~ /^PAD/) {
		$sample =~ /^(PAD\d+)(A|B|C)/; ### edit here ###
		$env = $1;
	}
	elsif ($sample =~ /^WC/) {
		$sample =~ /^(WC\d+)(A|B|C)/;
		$env = $1;
	}

	print OUT "$sample\t$env\n";

	$line=();
	@line=();
	$sample=();
	$env=();
}
close IN;
close OUT;

#!/usr/bin/perl
#June 6, 2012 by Terri Porter
#Script to parse a directory of fasta files to find number of baits/seqs in each
#usage perl fasta_stat_dir.plx

use strict;
use warnings;

#declare var
my $i=0;
my $file;
my $gb;
my $j=0;
my $line;
my $count=0;

#declare array
my @output;
my @in;
my @file;

@output = qx(ls | grep '.fastq\$\');

open (OUT, ">>", "stats.txt") || die "Error cannot open outfile: $!\n";

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	open (IN, "<", $file) || die "Error cannot open $file: $!\n";
	@in = <IN>;
	close IN;

	@file = split(/\./, $file);
	$gb = $file[0];
	
	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($line =~ /^@/) { ### edit here to search a fasta or fastq file! ###
			$count++;
		}
		$j++;
		$line=();
	}
	$j=0;
	
	print OUT "$gb\t$count\n";
	$count=0;
	$i++;
	$file=();
	@in=();
	@file=();
	$gb=();
}
$i=0;
close OUT;

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
my $flag=0;
my $cat="";
my $scalar;

#declare array
my @output;
my @in;
my @file;
my @cat;

@output = qx(ls | grep 'gappyout\$\');

open (OUT, ">>", "gappyout.stats") || die "Error cannot open outfile: $!\n";

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

		if ($line =~ /^>/) {
			if ($flag==0) {
				$count++;
				$flag=1;
			}
			elsif ($flag==1) {
				$count++;
				$flag=2;
			}
			elsif ($flag==2) {
				$count++;
			}
		}
		
		elsif ($flag==1) {
			$cat = $cat.$line;
			#print "$cat\n";
		}

		$j++;
		$line=();
	}
	$j=0;
	$flag=0;
	@cat = split(//, $cat);
	$scalar = scalar(@cat);
	#print "$scalar\n";
	
	print OUT "$gb\t$count\t$scalar\n";
	$count=0;
	$i++;
	$file=();
	@in=();
	@file=();
	$gb=();
	@cat=();
	$scalar=();
	$cat=();
}
$i=0;
close OUT;

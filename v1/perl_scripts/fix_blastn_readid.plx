#!/usr/bin/perl
#April 11, 2013 by Terri Porter
#Script to fix Nagissa's readid's in blastn file so that MEGAN imports it properly.  In the future, be sure to fix readids before the BLAST step.
#usage perl fix_blastn_readid.plx file.blastn

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $filename;

#declare array
my @in;


open (IN, "<", $ARGV[0]) || die "Error cannot open file.blastn: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0].".readidfixed";
open (OUT, ">>", $filename) || die "Error cannot open file.blastn.readidfixed: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^Query=/) {
		$line =~ s/length=\d+$//g;
		print OUT $line."\n";
	}
	else {
		print OUT $line."\n";
	}

	$i++;
	$line=();
}
$i=0;
close OUT;

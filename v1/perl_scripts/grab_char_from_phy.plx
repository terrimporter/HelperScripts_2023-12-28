#!/usr/bin/perl
#Sept. 20, 2012 by Terri Porter
#Script to grab all .phy files from directory, grab char from first line
#usage perl grab_char_from_phy.phx

use strict;
use warnings;

#declare var
my $dir;
my $i=0;
my $fileName;
my $filePath;
my $j;
my $line;
my $numChar;

#declare array
my @files;
my @in;
my @line;
my @numChar;

@files = qx(ls | grep .phy);

while ($files[$i]) {
	$fileName = $files[$i];
	chomp $fileName;
	#print "$fileName opened\n"; #test

	open (IN, "<", $fileName) || die "Error cannot open $fileName: $!\n";
	@in = <IN>;
	close IN;

	$line = $in[0]; #grab first line only
	chomp $line;
	@line = split(/\t/, $line);
	$numChar = $line[1]; #grab alignment length (number of characters) only
	#print "$numChar characters\n"; #test
	push(@numChar, $numChar);

	$line=();
	@line=();
	$numChar=();

	$i++;
	$fileName=();
	$filePath=();
	@in=();
}
$i=0;

open (OUT, ">", "numChar.txt") || die "Error cannot open outfile: $!\n";

while ($numChar[$i]) {
	$numChar=$numChar[$i];
	print OUT "$numChar\n";
	$i++;
	$numChar=();
}
$i=0;
close OUT;

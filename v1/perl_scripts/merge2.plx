#!/usr/bin/perl
#July 8, 2008 written to concatenate multiple .out files from the FastBlast.plx program and must be copied into the directory to be copied to work
#Nov.17, 2010 modify to work properly.
#Dec.13, 2010 modify to prompt for command line input
#usage $perl merge.plx

use strict;
use warnings;

#declare variables
my $dir;
my $file;
my $filename;
my $pathtofile;

#declare arrays
my @files;
my @new_array;

print "Enter full path to directory containing files to merge (with trailing '/' ):\n";
$dir = <STDIN>;
chomp $dir;

#$dir = "/home/terri/Saina/MID2_trimmed/agrep_out/MID2_TE_545R_parts/";
opendir DH, $dir;
@files = readdir (DH);

#open merged file in append mode
open (OUT, ">>merge.txt") ||die ("Cannot open new merged file");

foreach $file (@files) {
	$filename = $file;
	$pathtofile = $dir.$filename;
	#open the files one by one
	open (FH, $pathtofile) || die ("Cannot open file in array");
	#read the contents
	@new_array = <FH>;
	#store in the merged file
	print OUT @new_array;
	close (FH);
	}
close (OUT);


#!/usr/bin/perl
#March 19, 2013 by Terri Porter
#Script to grab BS# from missing.txt (from find_missing_bestTrees.plx) and find associated infile for RAxML and put into new directory called 'new'
#USAGE perl find_missing_infiles.plx missing.txt

use strict;
use warnings;
use File::Path;
use File::Copy;

#declare var
my $i=0;
my $filename;
my $output;
my $path = "new/";
my $dirs;
my $pattern;
my $to;

#declare array
my @in;
my @filename;

if (defined $ARGV[0] && length $ARGV[0]) {

	open (IN, "<", $ARGV[0]) || die "Error cannot open missing.txt:$!\n";
	@in = <IN>;
	close IN;
}
else {
	print "You must enter an argument like this:\n\tperl find_missing_bestTrees.plx missing.txt\n\n";
}

while ($in[$i]) {
	$filename = $in[$i];
	chomp $filename;
	@filename = split(/\./, $filename);
	$pattern = $filename[3];
#	print "pattern: $pattern\n";

	$output = qx(ls | grep $pattern);
	chomp $output;
#	print "output: $output\n";
	
	if (defined $output && length $output > 0) {

		if (! -d $path) {
			$dirs = eval { mkpath($path)};
			die "Failed to create $path: $@\n" unless $dirs;
		}
		$to = $path.$output;
		copy($output,$path) or die "Failed to copy $output: $!\n";
	}
	else {
		print "Couldn't find $filename\n";
	}

	$i++;
	$filename='';
	$output='';
}
$i=0;

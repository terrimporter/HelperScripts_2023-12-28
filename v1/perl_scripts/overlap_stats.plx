#!/usr/bin/perl
#October 4, 2016 by Terri Porter
#Script to calculate num, min, max, mean, mode for overlapping regions from SeqPrep -cd
#USAGE $perl overlap_stats.plx $infile

use strict;
use warnings;

#use lib '/home/tp45/modules';
use Statistics::Lite qw(:all);

#declare variables
my $dir;
my $file;
my $path;
my $i=0;
my $line;
my $length;
my $base;
my $num;
my $min;
my $max;
my $mean;
my $mode;

#declare array
my @files;
my @infile;
my @seq;
my @split;
my @length;

print "Enter name of directory with .paired.aln files including final '/': \n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Could not open '$dir': $!\n";

@files = readdir DIR;

print "Sample\tNumSeqs\tMin\tMax\tMean\tMode\n";

foreach $file (@files) {

	if ($file eq "." || $file eq "..") {
		next;
	}

	$path = $dir.$file;

	open (IN, "<", $path) || die "Error cannot read $file: $!\n";
	@infile = <IN>;
	close IN;

	while($infile[$i]){
		$line = $infile[$i];
		chomp $line; #error

		if ($line =~ /^MERG/) { #error
			$line =~ s/^MERG: //g; 
			push (@seq, $line);
		}
		else {
			$i++;
			next;
		}
		$i++;
	}
	$i=0;
	$line=();

	while ($seq[$i]){
		@split = split(//, $seq[$i]);
		$length = scalar(@split);
		push (@length, $length);
		@split =();#empty array
		$i++;
	}
	$i=0;

	$min = min (@length);
	$max = max (@length);
	$mean = mean (@length);
	$mode = mode (@length);
	$num = scalar(@seq);

	$base = substr $file,0,5;

	print $base."\t".$num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";
	@infile=();
	@seq=();
	@length=();
}

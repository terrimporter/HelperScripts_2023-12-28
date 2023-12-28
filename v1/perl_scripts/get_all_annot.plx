#!/usr/bin/perl
#March 30, 2011 by Terri Porter
#Script to parse clone_summary_best_rank.txt and grab lines that were annotated and passed % support filter (ex. > 95% NJ bootstrap support)
#Num lines now equals total number annotated at this rank, subtract total number of correct at this rank, difference is number incorrectly annotated with 95% support (Type I error)
#usage $perl get_all_annot.plx clone_summary_best_rank.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $filename;
my $outfilename;
my $support;

#declare arry
my @in;
my @line;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@in = <IN>;
close IN;

$filename = $ARGV[0];
$outfilename = $filename.".annot";

open (OUT,">>",$outfilename) || die ("Error cannot write to outfile: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /nil/) {
		$i++;
		next;
	}
	else {
		@line = split(/\t/,$line);
		$support = $line[2];
		
		if ($support >= 95) {
			print OUT "$line\n";
		}
	}
	$i++;
}

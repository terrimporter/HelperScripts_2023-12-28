#!/usr/bin/perl
#March 29,2011 by Terri Porter
#Script to remove duplicate lines
#usage $perl remove_duplicate_lines.plx infile

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $previous_line='nil';

#declare array
my @in;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@in = <IN>;
close IN;

open (OUT,">>", "dupliates_removed.txt") || die ("Error cannot write to outfile: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line eq $previous_line) {
		$i++;
		next;
	}	
	else {
		print OUT "$line\n";
		$previous_line = $line;
	}
	$i++;
}

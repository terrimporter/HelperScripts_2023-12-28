#!/usr/bin/perl
#March 31, 2011 by Terri Porter
#Script to get full gb from ITS.blast.parsed.all
#usage $perl get_gb_from_parsed_file.plx ITS.blast.parsed.all

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $hit_gb;

#declare array
my @in;

open (IN,"<",$ARGV[0]) || die ("Error cannot read from ITS.blast.parsed.all: $!\n");
@in = <IN>;
close IN;

open (OUT,">>","hit_gb.list") || die ("Error cannot write to hit_gb.list: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /\t\w+\|\w+\.\d+\|/) {
		$line =~ /\t\w+\|(\w+\.\d+)\|/;
		$hit_gb = $1;
		print OUT "$hit_gb\n";
	}
	else {
		print $line."\n";#troubleshoot pattern match
	}
	$i++;
}
close OUT;

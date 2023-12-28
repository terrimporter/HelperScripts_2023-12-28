#!/usr/bin/perl
#March 31, 2011 by Terri Porter
#Script to get full gb from ITS.blast.parsed.all
#usage $perl get_gb_from_parsed_file.plx ITS.blast.parsed.all
#modified to work with merged blastn file
#usage $perl get_gb_from_parsed_file2.plx merged.blastn

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $hit_gb;
my $flag=0;
my $query_gi;
my $query_gb;

#declare array
my @in;
my @line;

open (IN,"<",$ARGV[0]) || die ("Error cannot read from merged.blastn: $!\n");
@in = <IN>;
close IN;

open (OUT,">>","hit_gb.list") || die ("Error cannot write to hit_gb.list: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0) {
		if ($line =~ /Query=/) {
			$line =~ s/Query=//;
			$line =~ s/ //g;
			@line = split (/\|/,$line);
			$query_gi = $line[0];
			$query_gb = $line[1];
			print OUT "$query_gi|$query_gb\t";
			$flag=1;
		}
	}
	elsif ($flag==1) {
		if ($line =~ /^\w+\|\w+\.\d+\|\s+/) {
			$line =~ /^\w+\|(\w+\.\d+)\|/;
			$hit_gb = $1;
			print OUT "$hit_gb\n";
			$flag=0;
		}
	}
	$i++;
}
close OUT;

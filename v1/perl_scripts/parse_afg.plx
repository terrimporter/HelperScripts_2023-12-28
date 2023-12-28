#!/usr/bin/perl
#May 6,2011 by Terri Porter
#Script to parse velvet.afg file to get map for velvet_contig and velvet_readids
#usage parse_afg.plx velvet.afg

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $elements=0;
my $j="nil";
my $line2;
my $velvet_contig;
my $velvet_read;

#declare array
my @afg;
my @map;

open (IN,"<",$ARGV[0]) || die ("Error cannot read velvet.afg file: $!\n");
@afg = <IN>;
close IN;

open (OUT,">>","afg.parsed") || die ("Error cannot write to afg.parsed file :$!\n");

while ($afg[$i]) {
	$line = $afg[$i];
	chomp $line;
	if ($line =~ /{CTG/) {
		#print "found {CTG\n";
		$elements = scalar(@map);
		if ($elements >0) {
			print OUT "@map\n";
			@map=();
			$j = $i+1;
                        $line2 = $afg[$j];
                        chomp $line2;
                        if ($line2 =~ /iid:\d+/) {
                                $line2 =~ /iid:(\d+)/;
                                $velvet_contig = $1;
                                push (@map,$velvet_contig);
			}
		}
		else {
			$j = $i+1;
			$line2 = $afg[$j];
			chomp $line2;
			if ($line2 =~ /iid:\d+/) {
				$line2 =~ /iid:(\d+)/;
				$velvet_contig = $1;
				push (@map,$velvet_contig);
			}
		}
	}
	elsif ($line =~ /^\{TLE/) {
		$j = $i+1;
		$line2 = $afg[$j];
		chomp $line2;
		if ($line2 =~ /src:\d+/) {
			$line2 =~ /src:(\d+)/;
			$velvet_read = $1;
			push(@map,$velvet_read);
		}
	}
	$i++;
}
print OUT "@map\n"; #dont' forget to print last element in array
close OUT;

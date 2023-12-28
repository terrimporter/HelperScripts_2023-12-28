#!/usr/bin/perl
#Nov.16,2010 by Terri Porter
#Script to sort file.qual by the IDs in the sorted .fna and remove first 10 bases
#usage $perl sort_qual_by_agrepID.plx file.qual sorted_trimmed.fna

use strict;
use warnings;

#declare var
my $line;
my $id_to_match;
my $id;
my $x;
my $match=0;
my $i=1;
my $trimmed_phred;

#declare array
my @fna;
my @ids;
my @line;

open (QUAL,"<",$ARGV[0]) || die ("Error: $!\n");

open (FNA,"<",$ARGV[1]) || die ("Error: $!\n");
while (<FNA>) {
	$line = $_;
	chomp $line;
	if ($line =~ />/) {
		$line =~ />(\w{14})/;
		$id_to_match = $1;
		push(@ids, $id_to_match);
	}
}
close FNA;

open (OUT,">>","trimmed.qual") || die ("Error: $!\n");

while (<QUAL>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ /^>(\w{14})/;
		$id = $1;
		foreach $x (@ids) {
			if ($x eq $id) {
				$match=1;
			}
		}
	}
	elsif ($match==1) {
		@line = split(/ /,$line);
		while ($i<=10) {
			shift(@line);
			$i++;
		}
		$trimmed_phred = join(' ',@line);
		print OUT ">$id\n$trimmed_phred\n";
		$i=1;
		$match=0
	}
}
close QUAL;
close OUT;

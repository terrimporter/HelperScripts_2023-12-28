#!/usr/bin/perl
#July 8, 2011 by Terri Porter
#Script to sort ITS1 or ITS2 extracted sequences by length.  Keep only fragments greater than 100bp.
#usage perl filter_fasta_by_length.plx infile.fasta

use warnings;
use strict;

#declare var
my $i=0;
my $line;
my $header;
my $j;
my $next_line;
my $length;
my $sequence;
my $min_length;

#declare array
my @in;
my @next_line;

open (IN,"<",$ARGV[0]) || die ("Error reading infile: $!\n");
@in = <IN>;
close IN;

print "Please enter minimum sequence length (ex. 100):\n";
$min_length = <STDIN>;
chomp $min_length;

open (OUT,">>","ITS.filtered") || die ("Error creating outfile: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$j=$i+1;
		$next_line = $in[$j];
		chomp $next_line;
		@next_line = split(//,$next_line);
		$length = scalar(@next_line);

		if ($length >= $min_length) {
			$sequence = join('',@next_line);
			print OUT "$header\n$sequence\n";
			$sequence=();
		}
		$i+=2;
		@next_line=();
	}
	else {
		$i++;
	}
}

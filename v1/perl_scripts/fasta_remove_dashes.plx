#!/usr/bin/perl

# Jan. 10, 2019 by Terri M. Porter
# Script to remove dashes in an aligned FASTA file
# Usage perl fasta_remove_dashes.plx < file.fasta

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $filename;

# declare array
my @in;

# Read in fASTA file
open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@in = <IN>;
close IN;

# Create new outfile for reformatted FASTA
$filename = $ARGV[0].".out";

open (OUT, ">", $filename) || die "Cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		print OUT $line."\n";
		$i++;
		next;
	}
	else {
		$line =~ s/-//g; # remove dashes
		print OUT $line."\n";
		$i++;
		next;
	}
}
$i=0;
close OUT;

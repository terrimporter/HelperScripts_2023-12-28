#!/usr/bin/perl
#Sept. 20, 2011 by Terri Porter
#Script to check features.txt (ITS) and count only lines with features.  Report this number in paper as to the # annotated features retrieved from genbank.
#usage perl check_features.plx features.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;

#declare array
my @in;

open (IN,"<",$ARGV[0]) || die "Error cannot read infile: $!\n";
@in = <IN>;
close IN;

open (OUT,">>","features.txt.checked") || die "Error cannot write to outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^\s+/) {
		$i++;
		next;
	}
	elsif ($line =~ /^ID\tOrganism/) {
		$i++;
		next;
	}
	else {
		print OUT "$line\n";
	}
	$i++;
}
close OUT;

#then do $wc -l features.txt.checked

#!/usr/bin/perl
#Modified to remove genus and family references so can be trained to order rank
#March 21, 2012 by Terri Porter
#Script to remove genus info from FungiLSU_train_taxid.txt / mytaxa.txt so that trained rank is moved up to family
#usage perl removeGenusLine.plx infile.txt

#declare var
my $i=0;
my $line;

#declare array
my @in;

use warnings;
use strict;

open (IN,"<",$ARGV[0]) || die "Error reading infile.txt: $!\n";
@in = <IN>;
close IN;

open (OUT,">>","mytaxa_nogenus_nofamily.txt") || die "Error writing outfile.txt: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /\*genus/) {
		$i++;
		next;
	}
	elsif ($line =~ /\*family/) {
		$i++;
		next;
	}
	else {
		print OUT "$line\n";
	}
	$i++;
}
$i=0;

close OUT;


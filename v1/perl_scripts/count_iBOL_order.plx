#!/usr/bin/perl
#June 17, 2013 by Terri Porter
#Script to parse iBOL.txt to grab order, and count abundance of each in the 3.75-v1 data release package
#usage perl count_iBOL_order.plx iBOL.txt

use strict;
use warnings;

#declare line
my $line;
my $i=0;
my $order;
my $count;

#declare array
my @in;
my @line;

#declare hash
my %order; #indexed by order name

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

#parse iBOLD infile to grab and count order abundance
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if (length($line) > 0) {
		@line = split(/\t/,$line);
		$order = $line[4];

		if (exists $order{$order}) {
			$count = $order{$order};
			$count++;
			$order{$order} = $count;
		}
		else {
			$order{$order} = 1;
		}

	}
	$i++;
	$line=();
	@line=();
	$order=();
	$count=();
}
$i=0;

#print an outfile
open (OUT, ">>", "iBOL_order_ALL.txt") || die "Cannot open outfile: $!\n";

while ( ($order,$count) = each (%order) ) {
	print OUT "$order\t$count\n";
}
close OUT;


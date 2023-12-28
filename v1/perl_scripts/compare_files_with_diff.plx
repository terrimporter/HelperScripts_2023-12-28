#!/usr/bin/perl
#Aug. 14, 2013 by Terri Porter
#Script to compare all .sff files in a directory to each other, then use the unix diff command to ensure they are different
#usage perl compare_files.plx
#In the future, do this for all files from Shadi just in case

use strict;
use warnings;

#declare var
my $scalar;
my $compare;
my $output;
my $file;
my $length=0;
my $i=0;

#declare array
my @sff_files;
my @output;

#compare files, pairwise with no repetition, with unix diff command
@sff_files = qx(ls | grep sff);
$scalar = scalar(@sff_files);
#print "$scalar\n";
#print "@sff_files\n"; #WC1AX and WC8AZ are listed next to each other

open (OUT, ">>", "output.txt") || die "Error cannot open outfile: $!\n";

while ($scalar > 1) {
	$compare = shift(@sff_files);
	chomp $compare;
	$scalar = scalar(@sff_files);
#	print "$scalar\n";

	while ($sff_files[$i]) {
		$file = $sff_files[$i];
		chomp $file;

		@output = `diff "$compare" "$file"`; #backtickes same as qx() or system("");
		$output = join("", @output); #in case output is on more than one line
	 	$length = length($output);

		print OUT "$compare\t$file\t$output\n";

		if ($length == 0) {
			print "$compare and $file files do not differ\n";
		}

		$file=();
		$output=();
		@output=();
		$length=0;
		$i++;
	}
	$i=0;
	
	$compare=();
}
close OUT;

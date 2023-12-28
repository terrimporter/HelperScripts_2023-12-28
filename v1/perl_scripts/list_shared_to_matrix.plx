#!/usr/bin/perl
#Aug. 8, 2013 by Terri Porter
#Script to convert list.shared file from mothur to an abundance based and presence-absence matrix for [R] vegan package because number of columns exceeds that which can be processed by excel ~ 16K
#usage perl list_shared_to_matrix.plx file_list.shared

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $scalar;
my $maxindex;
my $newline;
my $j=0;
my $element;
my $newline2;

#declare array
my @in;
my @line;
my @trimmed;
my @newline;
my @transformed;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (ABUNDANCE, ">>", "abundance.csv") || die "Error cannot open outfile1: $!\n";

open (INCIDENCE, ">>", "incidence.csv") || die "Error cannot open outfile2: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^label/) {
		@line = split(/\t/,$line);
		$scalar = scalar(@line);
		$maxindex = $scalar-1;
		@trimmed = @line[3..$maxindex];
		$newline = join(",",@trimmed);
		print ABUNDANCE ",".$newline."\n";
		print INCIDENCE ",".$newline."\n";
	}

	else {
		@line = split(/\t/,$line);
		$scalar = scalar(@line);
		$maxindex = $scalar-1;
		@trimmed = @line[1,3..$maxindex];
#		print "@trimmed\n"; #test ok
		$newline = join(",",@trimmed);
		print ABUNDANCE "$newline\n";
#		print "$newline\n"; #test ok

		@newline = split(',',$newline);
#		print "@newline\n";#test ok
#		my $test = scalar(@newline);
#		print "number of elements: $scalar\n"; #test ok

		foreach $element (@newline) {

			if ($element =~ /^(PAD|WC)/) {
#				print "$element,";
				push(@transformed,$element);
			}

			elsif ($element =~ /^\d+/) {
				if ($element == 0) {
					push(@transformed,$element);
#					print "0,";
				}
				elsif ($element >= 1) {
					$element = 1;
					push(@transformed,$element);
#					print "1,";
				}
			}

			$j++;
			$element=();
		}
		$j=0;
#		print "\n";

		$newline2 = join(',',@transformed);
		print INCIDENCE "$newline2\n";
		@transformed=();

	}

	$i++;
	$line=();
	@line=();
	$scalar=();
	$maxindex=();
	@trimmed=();
	@newline=();
	$newline=();
	$newline2=();

}
$i=0;

close ABUNDANCE;
close INCIDENCE;

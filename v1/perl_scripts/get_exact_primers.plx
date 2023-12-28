#!/usr/bin/perl
#Nov.17, 2010 by Terri Porter
#Script to get the list of exact primers for a given degenerate primer sequence that contain ambiguities and/or Inosine
#Usage $perl get_exact_primers.plx < primer.txt > exact_primer.txt

use strict;
use warnings;

#declare variable
my $line;
my $name;
my $degenerate_primer;
my $i=0;
my $degenerate_base;
my $new_exact_primer;
my $old_exact_primer;
my $x;
my $exact_base;
my $y;
my $z;
my $scalar;

#declare array
my @line;
my @names;
my @degenerate_primers;
my @degenerate_bases;
my @new_exact_primers;
my @old_exact_primers;

while (<>) {
	$line = $_;
	chomp $line;
	@line = split(/\t/,$line);
	$name = $line[0];
	$degenerate_primer = $line[1];
	push(@names,$name);
	push(@degenerate_primers,$degenerate_primer);
}

while ($degenerate_primers[$i]) {
	$degenerate_primer = $degenerate_primers[$i];
	@degenerate_bases = split(//,$degenerate_primer);
print "{record start\n\tOriginal (degenerate) primer sequence: @degenerate_bases\n"; #test
my $test = scalar (@degenerate_bases);
print "\tPrimer length: $test\n"; #test

	foreach $degenerate_base (@degenerate_bases) {
		if ($degenerate_base =~ /(A|C|T|G)/i) {
			$exact_base = $degenerate_base;
			$scalar = scalar(@new_exact_primers);
			if ($scalar >= 1) {
				foreach $x (@new_exact_primers) {
					push(@old_exact_primers,$x);
				}
				@new_exact_primers = ();
				foreach $old_exact_primer (@old_exact_primers) {
					$new_exact_primer = $old_exact_primer.$exact_base;
					push(@new_exact_primers,$new_exact_primer);
				}
			}
			elsif ($scalar == 0){
				$new_exact_primer = $exact_base;
				push(@new_exact_primers,$new_exact_primer);
			}
			@old_exact_primers = ();
		}
		else {
			foreach $x (@new_exact_primers) {
				push(@old_exact_primers,$x);
			}
			@new_exact_primers=();
			foreach $old_exact_primer (@old_exact_primers) {
				if ($degenerate_base =~ /M/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /R/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /W/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /S/i) {
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base=~ /Y/i) {
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base=~ /K/i) {
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /V/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /H/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /D/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers,$new_exact_primer);
				}
				elsif ($degenerate_base =~ /B/i) {
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers, $new_exact_primer);
				}
				elsif ($degenerate_base =~ /(X|N|I)/i) {
					$new_exact_primer = $old_exact_primer."A";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."C";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."G";
					push(@new_exact_primers,$new_exact_primer);
					$new_exact_primer = $old_exact_primer."T";
					push(@new_exact_primers,$new_exact_primer);
				}
				else {
					print "Unknown nucleotide ambiguity character or other punctuation present in primer sequence.\n";
				}
			}
			@old_exact_primers = ();
		}
	}
	$y = $names[$i];
	$z = scalar(@new_exact_primers);
	print "\t$y ($z exact primers):\n";
	foreach $z (@new_exact_primers) {
		print "\t\t$z\n";
	}
	print "record end}\n";
	$i++;
	@new_exact_primers = ();
}

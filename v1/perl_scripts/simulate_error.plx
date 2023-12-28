#!/usr/bin/perl
#May 9, 2011 by Terri Porter
#Script to simulate pyroseuqencing error by randomly inserting a different base ever 1/100 bases for 1% error rate or 1/1000 bases for 0.1% error rate
#usage $perl simulate_error.plx file.fasta

use strict;
use warnings;

#declare var
my $random_number;
my $i=0;
my $line;
my $j=0;
my $base;
my $random_base;
my $return_probability;
my $return_base;

#declare array
my @in;
my @line;

open (IN,"<",$ARGV[0]) || die ("Error cannot read fasta.file:$!\n");
@in = <IN>;
close IN;

open (OUT,">>","error.fasta") || die ("Error cannot write error.fasta: $!\n");

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	if ($line =~ /^>/) {
		print OUT "$line\n";
	}
	else {
		@line = split(//,$line);

		while ($line[$j]) {
			$base = $line[$j];
			error_probability();			
			if ($return_probability==0) {
				print OUT $base;
			}
			else {
				if ($base eq "A") {
					insert_new_base_not_A();
					print OUT $return_base;
				}
				elsif ($base eq "C") {
					insert_new_base_not_C();
					print OUT $return_base;
				}
				elsif ($base eq "T") {
					insert_new_base_not_T();
					print OUT $return_base;
				}
				elsif ($base eq "G") {
					insert_new_base_not_G();
					print OUT $return_base;
				}
				else { # if it's N, just leave it as is
					print OUT $base;
				}
			}
			$j++;
		}
		$j=0;
		print OUT "\n";
	}
	$i++;
}
close OUT;

####################

sub error_probability {
	$random_number = rand; # random number from 0 to 1
	if ($random_number<=0.0001) { #0.01==1/100==1% error rate OR 0.0001==1/1000==0.01% error rate
		$return_probability = 1;
	}
	else {
		$return_probability = 0;
	}
}

####################

sub insert_new_base_not_A {
	$random_base = rand;
	if ($random_base <=0.33) {
		$return_base =  "C";
	}
	elsif ($random_base<=0.66) {
		$return_base = "T";
	}
	else {
		$return_base = "G";
	}
}

####################

sub insert_new_base_not_C {
	$random_base = rand;
	if ($random_base <= 0.33) {
		$return_base = "A";
	}
	elsif ($random_base <= 0.66) {
	       $return_base = "T";
      	}
	else {
	$return_base = "G";
	}
}

####################

sub insert_new_base_not_T {
	$random_base = rand;
	if ($random_base <= 0.33) {
		$return_base = "A";
	}
	elsif ($random_base <= 0.66) {
		$return_base = "C";
	}
	else {
		$return_base = "G";
	}
}

####################

sub insert_new_base_not_G {
	$random_base = rand;
	if ($random_base <= 0.33) {
		$return_base = "A";
	}
	elsif ($random_base <= 0.66) {
		$return_base = "C";
	}
	else {
		$return_base = "T";
	}
}


#!/usr/bin/perl
# Teresita M. Porter, Feb. 26, 2019

# Script to salvage R1.fastq.gz and R2.fastq.gz files with a different number of entries
# USAGE perl match_paired_end_raw_reads_fastq.plx filename_R1.fastq.gz filename_R2.fastq.gz

use strict;
use warnings;

# declare var
my $r1 = $ARGV[0];
my $r2 = $ARGV[1];
my $i=0;
my $j;
my $k;
my $l;
my $line;
my $id;
my $r1out = $ARGV[0].".NEW";
my $r2out = $ARGV[1].".NEW";
my $value;
my $header;
my $seq;
my $sep;
my $qual;
my $length_seq;
my $length_qual;

# declare array
my @R1;
my @R2;
my @line;

# declare hash
my %id; # key = id, value = 0

# Read in the compressed files
@R1 = `zcat $ARGV[0]`;
@R2 = `zcat $ARGV[1]`;

# Check number of lines in decompressed file
my $length_R1 = scalar @R1;
my $length_R2 = scalar @R2;
print "The R1 file has $length_R1 lines\nThe R2 file has $length_R2 lines\n\n";

# Hash the unique id from the fastq headers in the R1 file
while ($R1[$i]) {
	$line = $R1[$i];
	chomp $line;
	print $line."\n"; # test

	if ($line =~ /^\@/) {
#	print $line."\n"; # test
		@line = split(/ /,$line);
		$id = $line[0];

		# add unique ids to hash 
		$id{$id} = 0; # hashed, but no match with R2 made yet
		$i++;
		next;
	}
	else {
		$i++;
	}

}
$i=0;

# test check length of hash
my $size = keys %id;
print "The hash has $size entries\n";

# Check the unique id from the fastq headers in the R2 file
while ($R2[$i]) {
	$line = $R2[$i];
	chomp $line;

	if ($line =~ /^\@/) {
		@line = split(/ /,$line);
		$id = $line[0];

		if (exists $id{$id}) {
			$id{$id} = 1; # indicates matching id in R1 and R2
			$i++;
			next;
		}
		else {
			$i++;
			next;
		}
	}
	else {
		$i++;
		next;
	}
}
$i=0;

# Test check number of hash keys with R1 and R2 matches that are left now
while ( my ($key, $value) = each(%id)) {
	if ($value == "0") {
		delete $id{$key};
	}
}
$size = keys %id;
print "The hash has $size entries left after removing ones that don't match\n";

# Create a new R1 file

open (R1OUT, ">>", $r1out) || die "Error cannot open new R1 outfile: $!\n";
open (R2OUT, ">>", $r2out) || die "Error cannot open new R2 outfile: $!\n";

while ($R1[$i]) {
	$line = $R1[$i];
	chomp $line;

	if ($line =~ /\^@/ ) {
		@line = split(/ /,$line);
		$id = $line[0];

		if (exists $id{$id}) {
			$value = $id{$id};
			
			if ($value == 1) {
				$header = $line;

				$j = $i+1;
				$seq = $line[$j];
				chomp $seq;
				$length_seq = length $seq;
				if ($seq !~ /^[A-Z]/) { # seq must begin with A-Z letter or indicates messed up entry
					$i++;
					next;
				}
				$k = $j+1;

				$sep = $line[$k];
				chomp $sep;
				if ($seq !~ /^\+/) { # separator must be a plus sign or indicates messed up entry
					$i++;
					next;
				}
				$l = $k+1;

				$qual = $line[$l];
				chomp $qual;
				$length_qual = length $qual;
				if ($qual =~ /^\@/ || $qual =~ /^\+/) { # qual line must not start with @ or +
					$i++;
					next;
				}
				elsif ($length_seq != $length_qual) { # ensure length of sequence and quality scores is equal
					$i++;
					next;
				}
				else {
					print R1OUT $header."\n".$seq."\n".$sep."\n".$qual."\n";
				}
			}
			else {
				$i++;
				next;
			}
		}
		else {
			$i++;
			next;
		}
	}
	else {
		$i++;
		next;
	}
}
$i=0;
close R1OUT;

while ($R2[$i]) {
	$line = $R2[$i];
	chomp $line;

	if ($line =~ /^\@/ ) {
		@line = split(/ /,$line);
		$id = $line[0];

		if (exists $id{$id}) {
			$value = $id{$id};
			
			if ($value == 1) {
				$header = $line;

				$j = $i+1;
				$seq = $line[$j];
				chomp $seq;
				$length_seq = length $seq;
				if ($seq !~ /^[A-Z]/) { # seq must begin with A-Z letter or indicates messed up entry
					$i++;
					next;
				}
				$k = $j+1;

				$sep = $line[$k];
				chomp $sep;
				if ($seq !~ /^\+/) { # separator must be a plus sign or indicates messed up entry
					$i++;
					next;
				}
				$l = $k+1;

				$qual = $line[$l];
				chomp $qual;
				$length_qual = length $qual;
				if ($qual =~ /^\@/ || $qual =~ /^\+/) { # qual line must not start with @ or +
					$i++;
					next;
				}
				if ($length_seq != $length_qual) { # ensure length of sequence and quality scores is equal
					$i++;
					next;
				}
				
				print R1OUT $header."\n".$seq."\n".$sep."\n".$qual."\n";

			}
			else {
				$i++;
				next;
			}
		}
		else {
			$i++;
			next;
		}
	}
	else {
		$i++;
		next;
	}	
}
$i=0;
close R2OUT;

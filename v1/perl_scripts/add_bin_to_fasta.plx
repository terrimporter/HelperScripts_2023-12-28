#!/usr/bin/perl
# Teresita M. Porter, July 27, 2021
# BEFORE parsing files from BOLD, fix line endings ctrl-v ctrl-m, fix file format ff=unix, fix file encoding encoding=utf-8
# Script to add BOLD bin id's to the FASTA files from BOLD
# Remove terminal dashes, replace internal dashes with Ns
# USAGE perl add_bin_to_fasta.plx Master.map Master.fasta

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $processid;
my $bin;
my $outfile = "Master_bin.fasta";
my $newline;
my $newline_bin;

# declare arrays
my @fasta;
my @map;
my @line;
my @newline;

# declare hashes
my %map; # key = processid, value = bin

# mapping file is in this form
# processid\tsampleid\tbin_uri
open (IN, "<", $ARGV[0]) || die "Cannot open mapping file: $!\n";
@map = <IN>;
close IN;

# hash the mapping file for easy lookups
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	# fields are space-delimited 
	@line = split(/\t/, $line);
	$processid = $line[0];
	$processid =~ s/ //g; #handle weird space
	$bin = $line[2];
	if (length $bin) {
		$bin =~ s/ //g; #handle weird space
		if ($bin =~ /^\S+/) {
			$map{$processid} = $bin;
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

#print Dumper(\%map);

# fasta files from BOLD are in this form
# >processid|species name|amplicon|accession/sampleid
open (IN2, "<", $ARGV[1]) || die "Cannot open FASTA file: $!\n";
@fasta = <IN2>;
close IN2;

open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) { # FASTA header
		$newline = $line;
		$newline =~ s/^>//g;
		@newline = split(/\|/, $newline);
		$processid = $newline[0];

		if (exists $map{$processid}) {
			$bin = $map{$processid};
			$newline_bin = join "|", $line, $bin;
			print OUT $newline_bin."\n";
		}
		else { # assumes strict fasta file
			print OUT $line."|\n";
		}
	}
	else { # assumes strict FASTA format
		$line =~ s/^-+//g;
		$line =~ s/-+$//g;
		$line =~ s/-/N/g;
		print OUT $line."\n";
	}
	$i++;
	$bin=();

}
$i=0;

close OUT;


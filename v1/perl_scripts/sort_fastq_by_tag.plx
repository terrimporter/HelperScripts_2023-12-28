#!/usr/bin/perl

# Script to sort a fastq.gz file according to 5' and 3' tags (5bp), exact matches only for now, into their own fastq.gz files
# USAGE perl sort_fasta_by_tag.plx file.fastq.gz tags.txt

use strict;
use warnings;

# scalars
my $line;
my $tag;
my $ftag;
my $rtag;
my $rev;
my $revcomp;
my $cattag;
my $outfile;
my $file;
my $infile;
my $path;
my $filename;
my $i=0;
my $header;
my $j;
my $k;
my $l;
my $seq;
my $line3;
my $qual;

# arrays
my @tags;
my @line;
my @fastq;

# hashes
my %tags; #key = concatenated fwd and rev tag seqs, value = 1

# read in tag combinations

open (IN, "<", $ARGV[1]) || die "Cannot open infile: $!\n";
@tags = <IN>;
close IN;

foreach $tag (@tags) {
	$line = $tag;
	chomp $line;
	@line = split(/\t/,$line);
	$ftag = $line[1];
	$rtag = $line[2];
	$cattag = $ftag.$rtag;

	# assign to forward and reverse tag hashes 
	$tags{$cattag} = 1;

	# create a bunch of outifles to populate
	$outfile = $cattag.".fastq";
	open (OUT, ">", $outfile) || die "Cannot create outfile: $!\n";
	close OUT;
}

# read in fastq.gz file and process
$path = `pwd`;
chomp $path;
$file = $ARGV[0];
chomp $file;
$infile = $path."/".$file;
open (IN2, "gunzip -c $infile |" ) || die "Cannot open pipe to infile: $!\n";
@fastq = <IN2>;
close IN2;

while ($fastq[$i]) {
	$line = $fastq[$i];
	chomp $line;

	if ($line =~ /^@/) { #new fastq entry id line
		$header = $line; # get fastq header line

		$j = $i+1;
		$seq = $fastq[$j]; # get fastq seq line
		chomp $seq;

		# get fwd and rev tags
		$ftag = substr($seq, 0, 5);
		$rtag = substr($seq, -5, 5);

		#reverse complement the rtag
		$rev = reverse $rtag;
		$rev =~ tr/ATGCatgc/TACGtacg/;
		$revcomp = $rev;

		$cattag = $ftag.$revcomp;

		# get fastq line 3
		$k = $i + 2;
		$line3 = $fastq[$k];
		chomp $line3;

		# get fastq qual
		$l = $i + 3;
		$qual = $fastq[$l];
		chomp $qual;

		if (exists $tags{$cattag}) {
			$filename = $cattag.".fastq";
			open (IN3, ">>", $filename) || die "Cannot open fastq file: $!\n";
			print IN3 "$header\n$seq\n$line3\n$qual\n";
			close IN3;
			$i+=4; # go to next fastq entry
		}
		else { # go to next fastq entry
			$i+=4;
		}
	}
	else {
		$i++;
	}

}
$i=0;

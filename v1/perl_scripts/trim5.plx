#!/usr/bin/perl
#Nov. 18, 2011 by Terri Poter
#Edit to remove first 5 bp (barcode) from .fna, also remove first 5 phred scores from .qual
#NEW USAGE $perl trim5.plx
#Nov.16, 2010 by Terri Porter
#Script to remove first 10 base pairs from fasta file sorted by MID (10bp)
#usage $perl trim10.plx < fasta.file > trimmed_fasta.file

use strict;
use warnings;

#declare var
my $dir;
my $file;
my $filename;
my $path_to_file;
my $outfilename;
my $i=0;
my $line;
my $k=0;
my $header;
my $seq;
my $j=1;
my $trimmed_seq;
my $phredline;
my $trimmed_phredline;

#declare array
my @fna;
my @qual;
my @file;
my @seq;
my @phredline;

#get name of directory that contains .fna and .qual files to be trimmed
print "Enter path to directory (including final /):\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR,$dir) || die "Error cannot open directory: $!\n";
while ($file = readdir(DIR)) {
	next if ($file =~ /^\./); #skip the . and .. files
	if ($file =~ /\.fasta$/) {
		push(@fna,$file);
	}
	elsif ($file =~ /\.qual$/) {
		push (@qual,$file);
	}
}
closedir DIR;

#process .fna files first
while ($fna[$i]) {
	$filename = $fna[$i];
	$path_to_file = $dir.$filename;
	
	open (IN, "<", $path_to_file) || die "Cannot open file: $!\n";
	@file = <IN>;
	close IN;

	$outfilename = $filename.".trimmed";
	open (OUT,">>", $outfilename) || die "Cannot write to outfile: $!\n";
	
	while ($file[$k]) {#parse through each line of the file
		$line = $file[$k];
		chomp $line;
		
		if ($line =~ /^>/) {
			$header = $line;
		}
		else {
			$seq = $line;
			@seq = split (//,$seq);
			
			while ($j <= 5) {##### edit barcode length to remove here#####
				shift(@seq);
				$j++;
			}
			$trimmed_seq = join('',@seq);
			print OUT "$header\n$trimmed_seq\n";	
			$j=1;
		}
		$k++;
	}
	$k=0;
	close OUT;
	$i++;
}
$i=0;

#process .qual files next
while ($qual[$i]) {
	$filename = $qual[$i];
	$path_to_file = $dir.$filename;
	
	open (IN, "<", $path_to_file) || die "Cannot open file: $!\n";
	@file = <IN>;
	close IN;

	$outfilename = $filename.".trimmed";
	open (OUT, ">>", $outfilename) || die "Cannot write to outfile: $!\n";

	while ($file[$k]) {
		$line = $file[$k];
		chomp $line;
		
		if ($line =~ /^>/) {
			$header = $line;
		}
		else {
			$phredline = $line;
			@phredline = split (/ /,$phredline);
			
			while ($j <= 5) {##### edit barcode length to remove here too#####
				shift(@phredline);
				$j++;
			}
			$trimmed_phredline = join(' ',@phredline);
			print OUT "$header\n$trimmed_phredline\n";	
			$j=1;
		}
		$k++;
	}
	$k=0;
	close OUT;
	$i++;
}
$i=0;

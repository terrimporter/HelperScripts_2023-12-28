#!/usr/bin/perl
#Feb. 6, 2017 by Terri Porter
#Script to reformat new fasta headers from usearch9 unoise2 to look like the original format
#USAGE perl reformat_usearch9_fastaheaders.plx infile

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $id;
my $size;
my $outfilename;

#declare array
my @fasta;
my @line;

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@fasta = <IN>;
close IN;

$outfilename = $ARGV[0].".renamed";

open (OUT, ">>", $outfilename) || die "Cannot open outfile: $!\n";

while ($fasta[$i]){
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g;
		@line = split(/;/,$line);
		$id = $line[1];
		$id =~ s/uniq=//;
		$size = $line[3];
		print OUT ">".$id.";".$size.";\n";
	}
	else {
		print OUT $line."\n";
	}
	$i++;
}
$i=0;
close OUT;

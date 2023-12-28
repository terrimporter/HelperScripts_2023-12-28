#!/usr/bin/perl
#Dec. 18, 2013 by Terri Porter
#Script to trim 5bp tag + 22bp FC fwd primer and 20bp FC rev primer + 5bp tag for Shadi's barcode illumina sequences
#NEW USAGE perl trim5p3p.plx file.fasta
#use with ls | grep fasta | parallel -j 10 "perl trim5p3p.plx {}"
#
#Nov.16, 2010 by Terri Porter
#Script to remove first 10 base pairs from fasta file sorted by MID (10bp)
#usage $perl trim10.plx < fasta.file > trimmed_fasta.file

use strict;
use warnings;

#declare var
my $line;
my $header;
my $seq;
my $i=0;
my $trimmed_seq;
my $trimmed_seq2;
my $scalar;
my $length;
my $filename;
my $prefix;
my $newfile;

#declare array
my @in;
my @line;
my @trimmed_seq;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0];
chomp $filename;
$filename =~ /(\d+\D+)\.fasta/;
$prefix = $1;
$newfile = $prefix.".fasta.trimmed";

open (OUT, ">>", $newfile) || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/ /_/; #substitute first space with underscore
	}
	else {
		$seq = $line;
		$trimmed_seq = substr $seq, 25; #removes 5' 5bp tag + 22bp FC fwd primer; or 5' 5bp tag + 20bp BR fwd primer
		@trimmed_seq = split (//,$trimmed_seq);
		$scalar = scalar(@trimmed_seq);
		$length = $scalar - 23;
		$trimmed_seq2 = substr $trimmed_seq, 0, $length; #removes 3' 20bp FC rev primer + 5bp tag; or 3' 18bp BR rev primer + 5bp tag
		print OUT "$header\n$trimmed_seq2\n";
		$header=();
	}

	$i++;
	$line=();
#	$header=();
	$seq=();
	$trimmed_seq=();
	@trimmed_seq=();
	$scalar=();
	$length=();
	$trimmed_seq2=();

}
$i=0;
close OUT;

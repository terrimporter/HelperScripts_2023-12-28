#!/usr/bin/perl
#Apr. 10/17 edit to process invert orfs first
#Terri Porter, Mar. 29, 2017
#Script to filter fasta files by a file of orfs that meet certain criteria (start with invert mit genome, min length 309nt for BR5) when using ORFfinder
#USAGE perl filter_by_orf.plx infile.fasta

use strict;
use warnings;

#declare var
my $orffile;
my $outfile;
my $i=0;
my $line;
my $header;
my $j;
my $nextline;

#declare array
my @fasta;
my @orf;
my @headers;

#declare hash

open (FASTA, "<", $ARGV[0]) || die "Error cannot open fasta file: $!\n";
@fasta = <FASTA>;
close FASTA;

$orffile = $ARGV[0].".invert";

open (ORF, "<", $orffile) || die "Error cannot open orf file: $!\n";
@orf = <ORF>;
close ORF;

$outfile = $ARGV[0].".invertfiltered";

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$header = $line;
#		push (@headers,$header);
		
		if (grep /$header/, @orf) {
			#print "found match\n";#test
			print OUT ">".$line."\n";
			$j=$i+1;
			$nextline = $fasta[$j];
			chomp $nextline;
			print OUT $nextline."\n";
		}
	}
	$i++;
	$j=();
}
$i=0;
close OUT;

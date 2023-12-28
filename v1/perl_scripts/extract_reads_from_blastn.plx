#!/usr/bin/perl
#April 12, 2013 by Terri Porter
#Script to extract blastn results for a set of readids
#Get list of readids from run_usearch_OTU_clustering.plx
#Use new blast file for MEGAN analyses
#usage perl extract_reads_from_blastn.plx otus_minsize2.fa file.blastn.readidfixed

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $flag=0;
my $readid;
my $filename;

#declare array
my @otu;
my @blast;

#declare hash
my %readid; #indexed by readid

open (OTU, "<", $ARGV[0]) || die "Error cannot open readid list: $!\n";
@otu = <OTU>;
close OTU;

open (BLAST, "<", $ARGV[1]) || die "Error cannot open blast file: $!\n";
@blast = <BLAST>;
close BLAST;

#grab OTU readids
while ($otu[$i]) {
	$line = $otu[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g;
		$line =~ s/length=\d+;size=\d+;//g;
#		print $line."\n";#test
		$readid{$line} = 1;
	}
	$i++;
	$line=();
}
$i=0;

#extract BLAST reports
$filename = $ARGV[1].".extracted";
open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";

while ($blast[$i]) {
	$line = $blast[$i];
	chomp $line;
	
	if ($line =~ /^Query=/) {
		$flag=1;
	}

	if ($flag == 0) {
		print OUT $line."\n";
	}
	else {
		if ($line =~ /^Query=/) {
			$line =~ /^Query=\s{1}(\w+)/;
			$readid = $1;
			if (exists $readid{$readid}) {
				print OUT $line."\n";
				$flag=0;
			}
		}
	}
	
	$i++;
	$line=();
	$readid=();

}
$i=0;
close OUT;

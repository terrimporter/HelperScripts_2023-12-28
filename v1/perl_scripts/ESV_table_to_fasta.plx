#!/usr/bin/perl
# Teresita M. Porter, Jan. 7, 2021
# Script to use an ESV.table (based on exactly mapping primer-trimmed reads to denoised-chimera-free ESVs)
# to re-create a fasta file in a format useful for QIIME & PICRUSt (ex. >sample_Zotu)
# USAGE perl ESV_table_to_fasta.plx ESV.table cat.denoised.nonchimeras

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $otuID;
my $sample;
my $count;
my $outfile = "seqs.fna";
my $fullseq='';
my $newheader;
my $j=0;

# declare array
my @table;
my @samples;
my @counts;
my @fasta;
my @newheaders;

# declare hashes
my %table; #outer key = sample, inner key = sample, value = read count

open (IN, "<", $ARGV[0]) || die "Cannot open ESV.table :$!\n";
@table = <IN>;
close IN;

# hash ESV.table for easy lookup
while ($table[$i]) {
	$line = $table[$i];
	chomp $line;

	if ($line =~ /^#/) {# found headers
		@samples = split(/\t/,$line);
		$otuID = shift(@samples);
	}
	else {# found count values
		@counts = split(/\t/,$line);
		$otuID = shift(@counts);

		foreach $sample (@samples) {
			$count = shift(@counts);
			$sample =~ s/_//g;
			$table{$otuID}{$sample} = $count;
		}

	}
	$i++;
}
$i=0;


open (IN2, "<", $ARGV[1]) || die "Cannot open cat.denoised.nonchimeras: $!\n";
@fasta = <IN2>;
close IN2;

open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) { # found header
		if ($fullseq ne '') {
			if (length @newheaders > 0 ) {
				foreach $newheader (@newheaders) {
					print OUT $newheader."\n";
					print OUT $fullseq."\n";
				}
			}
			$fullseq='';
			$newheader='';
			@newheaders=();
		}
		
		$otuID = $line;
		$otuID =~ s/^>//;

		foreach $sample (@samples) {
			if (exists $table{$otuID}{$sample}) {
				$count = $table{$otuID}{$sample};
#				$sample =~ s/_//g; 
				$newheader = ">".$sample."_".$otuID;
				if ($count > 0) {
					while ($j < $count) {
						push(@newheaders,$newheader);
						$j++;
					}
					$j=0;
				}
			}
			else {
				print "Cannot find otuID $otuID sample $sample in hash\n";
			}
		}

	}

	else {
		$fullseq = $fullseq.$line;
	}

$i++;

}
$i=0;

# print last set of seqs
foreach $newheader (@newheaders) {
	print OUT $newheader."\n";
	print OUT $fullseq."\n";
}

close OUT;







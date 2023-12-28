#!/usr/bin/perl
#Tweak to convert USEARCH FASTA with line-breaks in sequence to a strict FASTA format
#Written by Terri Porter, lastest updat April 2, 2013, for Porter et al., 2013 MER
#Script to reformat fasta header to keep readid only
#USAGE $perl reformat_fsata_id.plx < infile > outfile

use strict;
use warnings;

#var
my $i=0;
my $line;
my $flag=0;
my $seq;
my $newseq;

#arrays
my @in;

@in = <STDIN>;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0 && $line =~ /^>/) { #first line in file
		$flag=1;
		print STDOUT $line."\n";
	}
	elsif ($flag==1) { #start the seq
		$seq = $line;
		$flag++;
	}
	elsif ($flag>0 && $line=~/^>/) { #print last sequence before going to next FASTA entry
		print STDOUT $seq."\n";
		$flag=1;
		print STDOUT $line."\n";
	}
	else { #build up the seq
		$newseq = $seq.$line;
		$seq = $newseq;
		$flag++;
	}
	$i++;
	$line=();
}
#don't forget to print the last seq
print STDOUT $seq."\n";
$i=0;

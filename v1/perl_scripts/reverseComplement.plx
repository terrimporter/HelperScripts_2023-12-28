#!/usr/bin/perl
#May 8, 2013 by Terri Porter
#Script to reverse-complement sequences in a fasta file
#usage perl reverseComplement.plx file.fasta
#OR use with parallel:
#ls | grep .fa | parallel -j 1 | 'perl reverseComplement.plx {}'

use strict;
use warnings;
use Bio::Perl;

#declare var
my $i=0;
my $flag=0;
my $line;
my $header;
my $j;
my $nextline;
my $originalSeq;
my $newSeq;
my $seq;
my $filename;
my $rc;
my $obj;

#declare array
my @fasta;

#declare hash
my %seq; #indexed by fasta header

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@fasta = <IN>;
close IN;

while($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if (length($line) > 0) { #skip empty lines

		if ($flag == 0 && $line =~ /^>/) {
			$header = $line;
			$j= $i+1;
			$nextline = $fasta[$j];
			chomp $nextline;

			if (length($nextline) > 0 ) {
			
				if ($nextline !~ /^>/) {
					$flag=1;
					$i++;
					next;
				}
				elsif ($nextline =~ /^>/) {
					$header=();
					$seq=();
					$originalSeq=();
					$newSeq=();
					$seq=();

					$flag=0;
					$i++;
					next;
				}
			}
			else {
				$i++;
				next;
			}
		}

		if ($flag == 1 && $line !~ /^>/) {

			if (exists $seq{$header}) {
				$originalSeq = $seq{$header};
				$newSeq = $originalSeq.$line;
				$seq{$header} = $newSeq;
			}
			else {
				$seq = $line;
				$seq{$header} = $seq;
			}
			$j = $i+1;
			$nextline = $fasta[$j];
			chomp $nextline;

			if (length($nextline) > 0) {
		
				if ($nextline !~ /^>/) {
					$flag=1;
					$i++;
					next;
				}
				elsif ($nextline =~ /^>/) {
					$header=();
					$originalSeq=();
					$newSeq=();
					$seq=();

					$flag=0;
					$i++;
					next;
				}
			}
			else {
				$i++;
				next;
			}
		}
	}
	else {
		$i++;
		next;
	}
}
$i=0;

$filename = $ARGV[0];
$filename = $filename.".rc";
open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";

while ( ($header,$seq) = each (%seq) ) {
	$rc = reverse_complement_as_string($seq);
	print OUT "$header\n$rc\n";
}
close OUT;

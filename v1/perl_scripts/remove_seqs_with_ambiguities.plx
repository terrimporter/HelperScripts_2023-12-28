#!/usr/bin/perl
# Teresita M. Porter, Feb. 1, 2020
# Read in strict FASTA file, remove sequences with any non-nucleotide characters (N's or ambiguities)
# USAGE perl remove_seqs_with_ambiguities.plx CO1v4BOLD.arth.fasta

use strict;
use warnings;

# vars
my $outfile = "CO1v4BOLD.arth.fasta.noambig";
my $outfile2 = "short550barcodes.fasta";
my $outfile3 = "ave550750barcodes.fasta";
my $outfile4 = "long750barcodes.fasta";
my $i=0;
my $line;
my $header;
my $j;
my $seq;
my $length;

# arrays
my @in; # FASTA file
my @seq;

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";
open (OUT2, ">>", $outfile2) || die "Cannot open outfile2: $!\n";
open (OUT3, ">>", $outfile3) || die "Cannot open outfile3: $!\n";
open (OUT4, ">>", $outfile4) || die "Cannot open outfile4: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) { # header
		$header = $line;
		$j = $i+1;
		$seq = $in[$j]; # seq
		chomp $seq;
		@seq = split(//,$seq);
		$length = scalar(@seq);
#		if ($length >= 600) { # 
			if ($seq !~ /[^ACGTacgt]/) { # screen out seqs with any non-nucleotide characters
				print OUT "$header\n$seq\n";

				if ($length < 550) {
					print OUT2 "$header\n$seq\n";
				}
				elsif ($length >=550 && $length <= 750) {
					print OUT3 "$header\n$seq\n";
				}
				elsif ($length > 750) {
					print OUT4 "$header\n$seq\n";
				}
			}
			else {
				$i++;
				next;
			}
#		}
#		else {
#			$i++;
#			next;
#		}
	}
	$i++;
}
$i=0;

close OUT;
close OUT2;
close OUT3;
close OUT4;

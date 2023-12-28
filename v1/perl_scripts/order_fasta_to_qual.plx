#!/usr/bin/perl
#Nov.18,2010 by Terri Porter
#Script to reorder nonredundant.fasta to match .qual
#usage $perl order_fasta_by_qual.plx file.fasta file.qual

use strict;
use warnings;

#declar var
my $i=0;
my $line;
my $id;
my $j=0;
my $line2;
my $id2;
my $k;
my $seq;

#declare array
my @fasta;
my @qual;


open (FASTA,"<",$ARGV[0]) || die ("Error: $!\n");
@fasta = <FASTA>;
close FASTA;

open (QUAL, "<",$ARGV[1]) || die ("Error: $!\n");
@qual = <QUAL>;
close QUAL;

open (OUT,">>","sorted.fasta") || die ("Error: $!\n");

while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ /(\w{14})/;
		$id = $1;
		$j=0;
		my $flag=0;
		while ($fasta[$j]) {
			$line2 = $fasta[$j];
			chomp $line2;
			if ($line2 =~ /^>/) {
				$line2 =~ /(\w{14})/;
				$id2 = $1;
				if ($id2 eq $id && $flag==0) {
					$k = $j+1;
					$seq = $fasta[$k];
					chomp $seq;
					print OUT ">$id2\n$seq\n";
					$flag=1;
				}
				if ($id2 eq $id && $flag==1){
					$j++;
					next;
				}
			}
			$j++;
		}
	}
	$i++;
}

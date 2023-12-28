#!/usr/bin/perl
#Dec. 5, 2011 by Terri Porter
#Script to consolidate ITS.filtered.reformatted from reformat_cluster_header3.plx
#with seqtrim.qual.newlineremoved.trimmed from trim_qual_by_outdata.plx
#to work with seqtrim.pl
#usage perl consolidate.plx file.fasta file.qual

use strict;
use warnings;
use Statistics::Lite qw(:all);
#use Data::Dump qw[pp];

#declare var
my $i=0;
my $line;
my $header;
my $id;
my $j;
my $seq;
my $qual;
my $mean;
my $num_elements;
my $line1;
my $line2;
my $x;
my $num;
my $check=();
my $length;

#declare array
my @fasta;
my @qual;
my @sorted_length;
my @sorted_qual;
my @phred;
my @unsorted;
my @sorted;
my @x;

#declare hash
my %fasta;
my %qual;
my %mean;
my %num;
my %map;
my %map_sorted;

open (FASTA,"<",$ARGV[0]) || die "Error cannot read from file.fasta: $!\n";
@fasta=<FASTA>;
close FASTA;

open (QUAL,"<",$ARGV[1]) || die "Error cannot read from file.qual: $!\n";
@qual=<QUAL>;
close QUAL;

#put fasta file into hash indexed by id
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w+)/;
		$id = $1;
		$j=$i+1;
		$seq = $fasta[$j];
		chomp $seq;
		$fasta{$id} = $seq;
	}
	$i++;
}
$i=0;

#put qual file into hash indexed by id
while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w+)/;
		$id = $1;
		$j=$i+1;
		$qual = $qual[$j];
		chomp $qual;
		$qual{$id} = $qual;

		@phred = split(/\s+/,$qual);
		$num_elements = scalar (@phred);
		$mean = mean (@phred);
		$mean{$id} = $mean;
		$num{$id} = $num_elements;
	}
	$i++;
}
$i=0;

#create array that includes id,length,avephred
while (($id,$num) = each (%num) ) {
	$mean = $mean{$id};
	$line = $id.",".$num.",".$mean;
	push(@unsorted,$line);
}

#Schwartzian Transform to sort multiple columns
@sorted = map { $_->[0]}
sort { $b->[2] <=> $a->[2] || $b->[3] <=> $a->[3] } 
map { [$_, split(/,/)] }
@unsorted;

open (FASTA,">>","file.fasta") || die "Error cannot write to fasta.file $!\n";
open (QUAL,">>","file.qual") || die "Error cannot write to fasta.qual $!\n";

foreach $x (@sorted) {
		@x = split(/,/,$x);
		$id = $x[0];
		if ($fasta{$id}) {
			$length = $num{$id};
			if ($length >= 50) { ##### only keep ITS2 if longer than 50 bp #####
				$seq = $fasta{$id};
				print FASTA ">$id\n$seq\n";
				$qual = $qual{$id};
				print QUAL ">$id\n$qual\n";
			}
		}
		$length=();
}

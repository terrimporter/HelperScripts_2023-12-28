#!/usr/bin/perl
#May 11, 2011 by Terri Porter
#Script to change fasta header from >gb_length_length_length to >gi|gb
#usage $perl fix_header.plx ITS_query.fasta gb_query.list gi_query.list

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb_noversion;
my $j=0;
my $gb_map;
my $k;
my $gi_map;
my $sequence;
my $filename;
my $filename2;

#declare array
my @fasta;
my @gb_map;
my @line;
my @gi_map;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile $!\n");
@fasta = <IN>;
close IN;

open (IN2,"<",$ARGV[1]) || die ("Error cannot read infile $!\n");
@gb_map = <IN2>;
close IN2;

open (IN3,"<",$ARGV[2]) || die ("Error cannot read infile $!\n");
@gi_map = <IN3>;
close IN3;

$filename = $ARGV[0];
$filename2 = $filename.".mapped";

open (OUT,">>",$filename2) || die ("Error cannot write outfile $!\n");

while ($fasta[$i]){
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		@line = split(/_/,$line);
		$gb_noversion = $line[0];

		while ($gb_map[$j]) {
			$gb_map = $gb_map[$j];
			chomp $gb_map;

			if ($gb_noversion eq $gb_map) {
				$gi_map = $gi_map[$j];
				chomp $gi_map;
				print OUT ">$gi_map|$gb_map\n";
				$k=$i+1;
				$sequence = $fasta[$k];
				chomp $sequence;
				print OUT "$sequence\n"; 
			}
			$j++;
		}
		$j=0;
	}
	$i++;
}
close OUT;

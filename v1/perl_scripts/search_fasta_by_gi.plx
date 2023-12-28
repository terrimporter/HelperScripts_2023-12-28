#!/usr/bin/perl
#Nov.7, 2011 by Terri Porter
#Script to parse a fasta file >gi|accession using a list of gi numbers
#usage perl search_fasta_by_gi.plx gb.list gi.list complete_overlap.gi file.fasta

use strict;
use warnings;

#declare array
my @gi;
my @fasta;
my @overlap;
my @gb;

#declare var
my $line;
my $i=0;
my $j=0;
my $k;
my $gi;
my $gb;
my $seq;
my $infile;
my $outfile;

#declare hash
my %map;

open (GB,"<",$ARGV[0]) || die "Error cannot read gi.list: $!\n";
@gb = <GB>;
close GB;

open (GI,"<",$ARGV[1]) || die "Error cannot read file.fasta: $!\n";
@gi = <GI>;
close GI;

#create mapping hash
while ($gi[$i]) {
	$gi = $gi[$i];
	chomp $gi;
	$gb = $gb[$i];
	chomp $gb;
	$gb =~ s/\.\d+//;
	$map{$gi} = $gb;
	$i++;
}
$i=0;

#test mapping hash
#while(my($key,$value) = each (%map)) {
#	print "$key\t$value\n";
#}

open (OVERLAP,"<",$ARGV[2]) || die "Error cannot read complete_overlap.gi: $!\n";
@overlap=<OVERLAP>;
close OVERLAP;

open (FASTA,"<",$ARGV[3]) || die "Error cannot read file.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

#automatically name outfile
$infile = $ARGV[3];
$outfile = $infile.".filtered";

open (OUT,">>","$outfile") || die "Error cannot write to outfile: $!\n";

#search for the fastas with complete overlap from four primers
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line=~ /^>/) {
		
		while ($overlap[$j]) {
			$gi = $overlap[$j];
			chomp $gi;
			$gb = $map{$gi};
			if ($line =~ /$gb/) {
				$k=$i+1;
				$seq = $fasta[$k];
				chomp $seq;

				print OUT "$line\n$seq\n";
			}
			$j++;
			$gi=();
			$gb=();
			$k=();
		}
		$j=0;
	}
	$i++;
}
$i=0;
close OUT;
			

#!/usr/bin/perl
#Sept.27, 2011 by Terri Porter
#Script to swap gb in full.fasta to gi to use with run_blastn_test.plx
#usage perl change_gb_to_gi.plx gb.query gi.query full.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $j=0;
my $gb;
my $short_gb;
my $gi;
my $line;
my $old_gb;
my $new_gi;
my $k;
my $seq;

#declare array
my @gb;
my @gi;
my @fasta;

#declare hash
my %map;

open (GB,"<",$ARGV[0]) || die "Error cannot read gb.query: $!\n";
@gb = <GB>;
close GB;

open (GI,"<",$ARGV[1]) || die "Error cannot read gi.query: $!\n";
@gi = <GI>;
close GI;

open (FASTA,"<",$ARGV[2]) || die "Error cannot read full.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

while ($gb[$i]) {
	$gb = $gb[$i];
	chomp $gb;
	$gb=~/(\w+)\.\d+/;
	$short_gb = $1;
	$gi = $gi[$i];
	chomp $gi;
	$map{$short_gb} = $gi;
	$i++;
}

#while (my($key,$value) = each (%map)) {
#	print "$key\t$value\n";#test
#}


open (MAPPED,">>","mapped_full.fasta") || die "Error cannot write to mapped_full.fasta: $!\n";

while ($fasta[$j]) {
	$line = $fasta[$j];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w+)/;
		$old_gb = $1;
		$new_gi = $map{$old_gb};
		print MAPPED ">$new_gi\n";
		$k=$j+1;
		$seq = $fasta[$k];
		chomp $seq;
		print MAPPED "$seq\n";
		$j++;
	}
	$j++;
}
close MAPPED;

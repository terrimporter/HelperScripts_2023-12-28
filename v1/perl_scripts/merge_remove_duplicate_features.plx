#!/usr/bin/perl
#Sept. 19, 2011 by Terri Porter
#Script to merge features_#.txt files and remove the duplicates
#usage perl merge_remove_duplicate_features.plx

use strict;
use warnings;

#declare var
my $file;
my $i=0;
my $j=0;
my $id;
my $features;
my $k=0;
my $line;
my $l=0;
my $key;
my $value;

#declare array
my @file_list;
my @current_file;
my @line;
my @unique;

#declare hash
my %features;

@file_list = qx(ls | grep "features_");

open (OUT,">>","features_merged_dereplicated.txt") || die "Error cannot print to outfile:$!\n";

while ($file_list[$j]) {
	$file = $file_list[$j];
	chomp $file;

	open (IN,"<",$file) || die "Error cannot read $file: $!\n";
	@current_file = <IN>;
	close IN;

	while ($current_file[$k]) {
		$line = $current_file[$k];
		chomp $line;
	
		if ($line =~ /^ID\t/) {
			$k++;
			next;
		}
		elsif ($line =~ /^(\w+)\t/) {
			$line =~ /^(\w+)\t(.+)$/;
			$id = $1;
			$features = $2;
			$features{$id} = $features;
		}
		$k++;
	}
	$k=0;
	$j++;
}

@unique = keys %features;
while ($unique[$l]) {
	$key = $unique[$l];
	$value= $features{$key};
	print OUT "$key\t$value\n";
	$l++;
}

close OUT;

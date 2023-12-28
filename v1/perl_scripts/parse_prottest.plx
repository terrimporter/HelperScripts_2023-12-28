#!/usr/bin/perl
#Oct. 6, 2011 by Terri Porter
#script to parse prottest outfiles to get best model according to AIC
#usage perl parse_prottest.plx

use strict;
use warnings;

#declare var
my $i=0;
my $infile;
my $ortho_id;
my $j=0;
my $line;
my $model;

#declare array
my @infiles;
my @file_contents;

open (OUT,">>","id_model.map") || die "Error cannot write to id_model.map: $!\n";

@infiles = qx(ls | grep .out);

while ($infiles[$i]) {
	$infile = $infiles[$i];
	chomp $infile;
	$infile =~ /^(\d+)\.phy.out/;
	$ortho_id = $1;

	open (IN,"<",$infile) || die "Error cannot read $infile:$!\n";
	@file_contents = <IN>;
	close IN;

	while ($file_contents[$j]) {
		$line = $file_contents[$j];
		chomp $line;

		if ($line =~ /^Best model according to AIC:/) {
			$line =~ /^Best model according to AIC:\s+(\S+)/;
			$model = $1;
			print OUT "$ortho_id\t$model\n";
		}
		$j++;
	}
	$j=0;
	$i++;
}
close OUT;

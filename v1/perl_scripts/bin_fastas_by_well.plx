#!/usr/bin/perl
#Dec. 10, 2013 by Terri Porter
#Script to sort sequences belonging to a single well into a single file for SS
#USAGE perl bin_fastas_by_well.plx PlateXX_derep_minsize1.derep.renamed.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $well;
my $filename;
my $j;
my $nextline;

#declare array
my @in;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in=<IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ />(\d+\D+)\s{1}/;
		$well = $1;
		
		#check if file exists else create new file
		$filename = $well.".fasta";
	
		open (FH, ">>", $filename) || die "Couldn't open $filename: $!\n";
		print FH $line."\n";
		$j=$i+1;
		$nextline = $in[$j];
		chomp $nextline;
		print FH $nextline."\n";
		
	}
	$i++;
	$line=();
	$well=();
	$filename=();
	$j=();
	$nextline=();
}
$i=0;


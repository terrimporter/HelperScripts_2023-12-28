#!/usr/bin/perl
#April 4, 2017 by Terri Porter
#Script to add read abundance from cat.table to cat.denoised.out
#USAGE perl add_abundance_to_rdp_out.plx cat.table cat.denoised.out

use strict;
use warnings;

#declare var
my $i=0;
my $global_otu;
my $line;
my $otupart;
my $assignment;
my $otu;
my $j=0;
my $sample;
my $abund;
my $outfile = "cat.denoised.out.updated";
my $newline;

#declare array
my @table;
my @rdp;
my @line;
my @otupart;
my @headers;
my @otu;

#declare hash
my %assignment; #key = global OTU, value = rdp taxonomic assignment
my %table; #hash of hashes, key1=otu, key2=sample, value=readnumber

open (IN, "<", $ARGV[0]) || die "Error cannot open infile1: $!\n";
@table = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Error cannot open infile2: $!\n";
@rdp = <IN2>;
close IN2;

#parse through global rdp taxonomic assignment
while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$otupart = shift @line;
	@otupart = split(/;/, $otupart);
	$global_otu = shift @otupart;
	$assignment = join("\t", @line);
	$assignment{$global_otu} = $assignment;
	$i++;
	@line=();
	@otupart=();

}
$i=0;

#grab list of headers and otus and build hash of hashes
while ($table[$i]) {
	$line = $table[$i];
	chomp $line;

	if ($i==0) { #header row
		@headers = split(/\t/,$line);
		shift @headers;
	}
	else {
		@line = split(/\t/,$line);
		$global_otu = shift @line;

		if (scalar @headers == scalar @line) {
#			print "Check PASSED\n"; #TEST OK
		}
		else {
			print "Check FAILED\n";
		}
		
		while ($line[$j]) {
			$table{$global_otu}{$headers[$j]} = $line[$j];
			$j++;
		}
		$j=0;

	}
	@line=();
	$i++;
}
$i=0;

#loop through hash of hashes, if OTU abund >= 3 keep and append taxonomic assignment, print out new assignment report


open (OUT, ">>", $outfile) || die "Error cannot open outfile : $!\n";

#add read abundance to rdp out as a new column
foreach $global_otu (sort keys %table) {
	foreach $sample (keys %{$table{$global_otu}}) {
		$abund = $table{$global_otu}{$sample};
		
		if ($abund >= 3) {
			$assignment = $assignment{$global_otu};
			print OUT "$global_otu\t$sample\t$abund\t$assignment\n";
		}
		else {
			next;
		}
	}
}

close OUT;

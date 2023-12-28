#!/usr/bin/perl

#Dec. 15, 2016 by Terri Porter
#Script to use to compare original OTUheaders with those after running blasts (max_hsps 1, max_target_seqs 1) and indicating which OTUheaders are missing because there were no BLAST hits with an evalue better than '1e-20'
#USAGE perl flag_no_BLAST_hits matrix_sequences_renamed.fasta matrix_sequences_renamed.fasta.blastn.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $cluster;
my $header;
my $value;

#declare array
my @fasta;
my @blast;
my @line;
my @allMatches;

#declare hash
my %cluster; #key=$cluster, value="1";

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@fasta = <IN>;
close IN;

#populate hash with complete set of clusterIDs
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g;
		@line = split(/;/,$line);
		$cluster = $line[0];
		$cluster =~ s/Shadi_CO1_//g; #clean it up
		$cluster{$cluster} = "0";
	}
	else {
		$i++;
		next;
	}
	$i++;
}
$i=0;

open (IN2, "<", $ARGV[1]) || die "Error cannot open infile2: $!\n";
@blast = <IN2>;
close IN2;

open (OUT, ">>", "matrix.blastn.parsed") || die "Error cannot open outfile: $!\n";

while ($blast[$i]) {
	$line = $blast[$i];
	chomp $line;

	if ($line =~ /^cluster/) {
		@line = split(/\t/,$line);
#		$header = $line[0];
#		@header = split(/;/,$header);
		$cluster = $line[0];
#		$cluster =~ s/Shadi_CO1_//g;

		if (exists $cluster{$cluster}) {
			print OUT $line."\n";
			$cluster{$cluster} = "1";
		}
		else {
			print OUT $header."\tNO BLAST hits\n";
		}
	}	
	else {
		$i++;
		next;
	}
	$i++;

}
$i=0;

@allMatches = grep { $cluster{$_} eq "0" } keys %cluster;

foreach $cluster (@allMatches) {
	print OUT $cluster."\tNo BLAST matches\n";
}

close OUT;

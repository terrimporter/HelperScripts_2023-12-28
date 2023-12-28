#!/usr/bin/perl
#March 15,2011 by Terri Porter
#Script to use a mapped fasta file (with gi numbers in the header) to conduct sap searches using the forceexcludegilist option.
#usage $perl run_sap_test.plx ITS.fasta.mapped
#modified to run blastn (BLAST+) repeated using a 'leave one out' approach

#declare var
my $line;
my $i=0;
my $gi;
my $filename;
my $j=0;
my $outfilename;
my $gifilename;
#declare array
my @fasta;
my @line;
my @gi;
my @output;

use strict;
use warnings;

open (FASTA,"<",$ARGV[0]) || die ("Error cannot read from mapped fasta file: $!\n");
@fasta = <FASTA>;
close FASTA;

print "\nCreating individual gi.fasta files\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ s/^>//;
		@line = split(/\|/,$line);
		$gi = $line[0];
		push(@gi,$gi);
		$filename = $gi.".fasta";
		$gifilename = $gi.".list";
		
		open (OUT,">>",$filename) || die ("Error cannot write to fasta outfile: $!\n");
		print OUT ">$line\n";
		
		open (OUT2,">>",$gifilename) || die ("Error cannot write to gi list outfile: $!\n");
		print OUT2 "$gi\n";
		
	}
	else {
		print OUT "$line\n";
		close OUT;
		close OUT2;
	}
	$i++;
}

print "\nStarting blastn searches\n";

while ($gi[$j]) {
	$gi = $gi[$j];
	$filename = $gi.".fasta";
	$outfilename = $filename.".blastn";
	$gifilename = $gi.".list";
	my @output = qx(blastn -task blastn -db /net/bioinfo/1/blast/blastdb/nt -query $filename -out $outfilename -outfmt 0 -negative_gilist $gifilename);
	@output=();
	$j++;
}

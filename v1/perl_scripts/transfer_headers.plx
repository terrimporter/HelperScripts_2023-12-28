#!/usr/bin/perl
#Terri Porter, Jan. 18, 2013
#Script to copy headers from Ustilago maydis genes file and transfer to transcript file
#usage perl transfer_headers.plx genes.fasta transcripts.fasta
#be sure to remove newline from fasta part of files first!

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $header;
my $newline;

#declare array
my @genes;
my @transcripts;

open (GENES, "<", $ARGV[0]) || die "Error cannot read genes.fasta: $!\n";
@genes = <GENES>;
close GENES;

open (TRANSCRIPTS, "<", $ARGV[1]) || die "Error cannot read transcripts.fasta: $!\n";
@transcripts = <TRANSCRIPTS>;
close TRANSCRIPTS;

open (OUT, ">>", "new.fasta") || die "Error cannot read outfile: $!\n";

while ($genes[$i]) {
	$line = $genes[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$header = $1;
		print OUT ">$header\n";
	}
	else {
		$newline = $transcripts[$i];
		chomp $newline;
		print OUT "$newline\n";
	}
	$i++;
}
$i=0;

close OUT;


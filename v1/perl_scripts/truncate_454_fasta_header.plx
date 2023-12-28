#!/usr/bin/perl
#Terri Porter, Jan. 14, 2013
#Script to remove everything except 14 digit readid from 454 fasta files from Joel
#usage perl truncate_454_fasta_header.plx Joel.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
#my $part1;
my $readid;
my $filename;

#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open file.fasta: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0].".truncated";

open (OUT, ">>", $filename) || die "Error cannot open $filename: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/=/, $line);
		$readid = $line[0];
#		$readid =~ s/^>//g;
		$readid =~ s/rank//g;
		$readid =~ s/;//g; #remove any ; if present at the end of the line

		print OUT "$readid\n";
	}
	elsif ($line =~ /^\s+/) { #skip over any blank lines
		next;
	}
	else {
		print OUT "$line\n";
	}

	$i++;
	$line=();
	@line=();
#	$part1=();
	$readid=();
}
$i=0;



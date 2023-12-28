#!/usr/bin/perl
#Nov.23,2010 by Terri Porter
#Script to convert fasta file (from mothur) into a relaxed phylip format for RAxML
#usage $perl fasta_to_relaxedphylip.plx < mothur.fasta

use strict;
use warnings;

#declare var
my $line;
my $header;
my $id;
my $num_seq=0;
my $seq;
my $num_char;
my $x;
my $scalar;
my $sum=0;
my $average;
my $i=0;

#declare array
my @ids;
my @seqs;
my @split;
my @scalar;

while (<>) {
	$line = $_;
	chomp $line;
	if ($line = /^>/) {
		$header = $line;
		print $header."\n";
		$header =~ /\|(\w{14})/;
		$id = $1;
		push(@ids,$id);
		$num_seq++;
	}
	else {
		$seq = $line;
		print "$seq\n"; #test
		#$seq =~ s/\./N/g;
		push(@seqs,$seq);
		@split = split(//,$seq);
		$scalar = scalar(@split);
		push(@scalar,$scalar);
	}
}
print "@seqs\n"; #test

($sum += $_) for @scalar;

$average = $sum/$num_seq;

print "numid = $num_seq\nnumchar = $average\n";

open (OUT,">>","phylip.txt") || die ("Error: $!\n");

print OUT "$num_seq\t$average\n";
while ($ids[$i]) {
	$id = $_;
	$seq = $seqs[$i];
	print OUT "$id\t$seq\n";
	$i++;
}

#!/usr/bin/perl
#Oct.14,2010 by Terri Porter
#Script to sort seqtrim .fna and .qual by phred score prior to dereplication with Uclust
#usage $perl sort_seqtrim_by_qual.plx file.fna file.qual

use strict;
use warnings;

#declare variables
my $i=0;
my $line;
my $j=0;
my $phred_score;
my $k=0;
my $l=0;
my $index;
my $header;
my $phred_seq;
my $m=0;
my $fna_header;
my $fna_seq;
my $key;
my $value;

#declare arrays
my @fna;
my @qual;
my @header;
my @phred_seq;
my @phred_scores;
my @phred_count;
my @value;
my @fna_header;
my @fna_seq;

#declare hash
my %hash;

open (FNA,"<",$ARGV[0]) || die ("Error:$!\n");
@fna = <FNA>;
close FNA;

open (QUAL,"<",$ARGV[1]) || die ("Error:$!\n");
@qual = <QUAL>;
close QUAL;

while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		push(@header,$line);
	}
	else {
		push (@phred_seq,$line);
	}
	$i++;
}

#test
my $test = scalar(@header);
my $test2 = scalar(@phred_seq);
print "header: $test\nphred_seq: $test2\n";

while ($phred_seq[$j]) {
	$line = $phred_seq[$j];
	@phred_scores = split(/ /,$line);
	foreach $phred_score (@phred_scores) {
		if ($phred_score < 20) {
			$k++;
		}
	}
	#push(@phred_count,$a); #use hash instead to permit easier sorting and key/index retrieval
	$hash{$j} = $k;
	@phred_scores=();
	$k=0;
	$j++;
}

#sort by hash values
foreach $value (sort {$hash{$a} <=> $hash {$b}} keys %hash) {
	#print "$value $hash{$value}\n";
	push(@value,$value);
}

#test
my $test3 = scalar(@value);
print "value: $test3\n";

open (QUAL2,">>","sorted.qual") || die ("Error:$!\n");

foreach my $l (@value) {
	$index = $l;
	$header = $header[$index];
	$phred_seq = $phred_seq[$index];
	print QUAL2 "$header\n$phred_seq\n";
}
close QUAL2;

while ($fna[$m]) {
	$line = $fna[$m];
	chomp $line;
	if ($line =~ /^>/) {
		push(@fna_header,$line);
	}
	else {
		push(@fna_seq,$line);
	}
	$m++;
}

#test
my $test4 = scalar(@fna_header);
my $test5 = scalar(@fna_seq);
print "fna_header: $test4\nfna_seq: $test5\n";

open (FNA2,">>","sorted.fna") || die ("Error:$!\n");

foreach my $l (@value) {
	$index = $l;
	$fna_header = $fna_header[$index];
	$fna_seq = $fna_seq[$index];
	print FNA2 "$fna_header\n$fna_seq\n";
}
close FNA2;

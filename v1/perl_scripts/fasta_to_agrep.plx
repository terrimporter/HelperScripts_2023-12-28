#!/usr/bin/perl
#Nov.12,2010 by Terri Porter
#Script to reformat a fasta file so that it can be properly parsed by agrep.
#Use output to search for 5' primer sequences with/without mismatches allowed
#usage $perl fasta_to_agrep.plx file.fasta primer.txt
#Nov.24,2010 modify to include a remove_newline.plx part

use strict;
use warnings;

#declare var
my $line;
my $id;
my $seq;
my $primername;
my $primersequence;
my $i=0;
my $pattern;
my $fname;
my $filename;
my $k=0;

#declare array
my @output;
my @line;
my @primername;
my @primersequence;

open (IN1,"<",$ARGV[0]) || die "Error: $!\n";
open (TEMP1,">>","temp1.txt") || die ("Error:$!\n");

while (<IN1>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		print TEMP1 "\n",$line,"\n";
	}
	else {
		print TEMP1 $line;
	}
}

close IN1;
close TEMP1;

open (IN2,"<","temp1.txt") || die ("Error: $!\n");
open (TEMP2,">>","temp2.txt") || die ("Error:$!\n");

while (<IN2>) {
	$line=$_;
	chomp $line;
	if ($line =~ /\S+/) {
		print TEMP2 $line,"\n";
	}
}
close IN2;
close TEMP2;

open (IN3,"<","temp2.txt") || die ("Error: $!\n");
open (TEMP3,">>","temp3.txt") || die "Error: $!/n";

while (<IN3>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ /^>(\w{14})\s+/;
		$id = $1;
	}
	elsif ($line =~ /^\w+/) {
		$seq = $line;
		print TEMP3 "$line|$id\n";
	}
}
close IN3;
close TEMP3;

open (IN4,"<",$ARGV[1]) || die "Error: $!\n";

while (<IN4>) {
	$line = $_;
	chomp $line;
	@line = split(/\t/,$line);
	$primername = $line[0];
	$primersequence = $line[1];
	push(@primername,$primername);
	push(@primersequence,$primersequence);
}
close IN4;

while ($primersequence[$i]) {
	$pattern = $primersequence[$i];
	###adjust agrep settings here###
	@output = qx(agrep -i -1 -D1 -I1 "^$pattern" temp3.txt);
	$fname = $primername[$i];
	$filename = $fname.".agrep";
	open (OUT,">>",$filename) || die "Error: $!\n";
	while ($output[$k]) {
		$line = $output[$k];
		print OUT "$line";
		$k++;
	}
	close OUT;
	$k=0;
	$i++;
}

unlink("temp1.txt","temp2.txt","temp3.txt");

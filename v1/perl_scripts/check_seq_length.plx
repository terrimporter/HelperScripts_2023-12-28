#!/usr/bin/perl
#Nov. 22, 2016 by Terri Porter
#Script to verify that fragments are exactly 200bp long, if not, then do not use
#USAGE perl check_seq_length.plx primerName.fasta200

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $header;
my $length;
my $primerName;
my $outfile;

#declare array
my @in;
my @seq;
my @basename;

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@in = <IN>;
close IN;

@basename = split(/\./,$ARGV[0]);
$primerName = $basename[0];

#print "primer name: $primerName\n"; #test

$outfile = $primerName.".fasta200checked";

open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;

#print "header: $header\n"; #test

		$i++;
		next;
	}
	else {
		@seq = split("",$line);
		$length = scalar @seq;

#print "length: $length\n"; #test

		if ($length == 200 ) {
			print OUT $header."\n".$line."\n";
			$header=();
			@seq=();
			$i++;
			next;
		}
		else {
			$header=();
			@seq=();
			$i++;
			next;
		}
	}
}
$i=0;
close OUT;

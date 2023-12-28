#!/usr/bin/perl
#Sept.26, 2011 by Terri Porter
#Script to change primer_400.fasta into shorter fragments of 50, 100, or 200 bp
#usage perl truncate_seqs.plx primer_400.fasta

use strict;
use warnings;

#declare var
my $trim_from;
my $frag_size;
my $i=0;
my $line;
my $j;
my $header;
my $sequence;
my $truncated_seq;
my $original_length;
my $length_to_keep;
my $offset;
my $truncated_seq;

#declare array
my @in;
my @sequence;

open (IN,"<",$ARGV[0]) || die "Error cannot read primer_400.fasta: $!\n";
@in = <IN>;
close IN;

print "Trim from the 5' or 3' end?  (enter 5 or 3):\n";
$trim_from = <STDIN>;
chomp $trim_from;

print "Enter fragment size to trim to:\n";
$frag_size = <STDIN>;
chomp $frag_size;

open (OUT,">>","truncated.fasta") || die "Error cannot write to truncated.fasta: $!\n";

while($in[$i]) {
	$line = $in[$i];
	chomp $line;
	$j=$i+1;

	if ($line =~ /^>/) {
		$header = $line;
		$sequence = $in[$j];
		chomp $sequence;
		
		@sequence = split(//,$sequence);
		$original_length = scalar(@sequence);
		#$length_to_keep = $original_length-$frag_size;
		
		if ($trim_from == 5) {
			$offset = $original_length-$frag_size;
			$truncated_seq = substr $sequence,$offset,$frag_size;
		}

		elsif ($trim_from == 3) {
			$truncated_seq = substr $sequence,0,$frag_size;
		}
		print OUT "$header\n$truncated_seq\n";
		$i++;
	}
	$i++;
	$header=();
	$sequence=();
	@sequence=();
	$original_length=();
	$length_to_keep=();
	$truncated_seq=();
	$offset=();
}
close OUT;

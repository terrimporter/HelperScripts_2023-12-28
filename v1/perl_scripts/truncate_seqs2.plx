#!/usr/bin/perl
#Nov.22, 2016 edit to provide all var on the command line
#Sept.26, 2011 by Terri Porter
#Script to change primer_400.fasta into shorter fragments of 50, 100, or 200 bp
#usage perl truncate_seqs.plx primer_400.fasta
#NEW usage perl truncate_seqs2.plx infile.fasta 5 200

use strict;
use warnings;

#declare var
my $filename;
my $primerName;
my $ext;
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
#my $truncated_seq;
my $outfile;

#declare array
my @in;
my @sequence;

open (IN,"<",$ARGV[0]) || die "Error cannot read primer_400.fasta: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0];
chomp $filename;
($primerName,$ext) = split(/\./,$filename);

#print "Trim from the 5' or 3' end?  (enter 5 or 3):\n";
$trim_from = $ARGV[1];
chomp $trim_from;

#print "Enter fragment size to trim to:\n";
$frag_size = $ARGV[2];
chomp $frag_size;

$outfile = $primerName.".fasta200";

open (OUT,">>",$outfile) || die "Error cannot write to primerName.fasta200: $!\n";

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

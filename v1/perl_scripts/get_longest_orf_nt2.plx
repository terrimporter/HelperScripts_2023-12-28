#!/home/terri/miniconda3/envs/myenv/bin/perl
# Sept 5., 2019 by Teresita M. Porter
# Script to grab longest CDS from ORFfinder nucleotide output (option 1)
## NOT using cutoff = minimum number of base pairs required to retain CDS ##
# print to STDOUT here, in snakemake, redirect to an outfile
# USAGE perl get_longest_orf.plx orf_nt.fasta

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $orf;
my $otu;
my $length;
my $flag=0;
my $seq;
my $longest;
#my $minimum = $ARGV[1]; #mode length after primer trimming BR5 = 310bpp use 309, F230R = 229bp use 219
#my $filename = "longest_".$minimum."_nt.fasta";

# declare array
my @in;
my @line;
my @otus;
my @longest;

# declare hash
my %length;
my %seq;

open (IN, "<", $ARGV[0]) || die "Can't open infile: $!\n";
@in=<IN>;
close IN;

# parse file to get fasta header stats first
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0 && $line =~/^>/) {
		($orf, $otu, $length) = parse_header($line);
		$length{$otu}{$orf} = $length; #create hash of hashes
		$flag = 1;
		$i++;
		next;
	}
	elsif ($flag==1) {
		$seq = $line;
		$flag = 2;
		$i++;
		next;
	}
	elsif ($flag==2 && $line !~ /^>/) {
		$seq = $seq.$line;
		$i++;
		next;
	}
	else {
		$flag = 1;
		$seq{$otu}{$orf} = $seq;
		($orf, $otu, $length) = parse_header($line);
		$length{$otu}{$orf} = $length;
		$i++;
		next;
	}
}
$i=0;
		
# add last seq to hash
$seq{$otu}{$orf} = $seq;

# loop through each otu, get orf id for longest length, print otu and seq to outfile

# get unique otu keys
@otus = keys %length;



# print longest ORF without doing any size filtering yet
while ($otus[$i]) {
	$otu = $otus[$i];

	@longest =  sort { $length{$otu}{$a} <=> $length{$otu}{$b} } keys %{$length{$otu}};
	$longest = $longest[-1]; #inner key to longest is at bottom of array
	$length = $length{$otu}{$longest};
#	if ($length >= $minimum) {
		$seq = $seq{$otu}{$longest};
		print STDOUT ">$otu\n$seq\n";
		$i++;
		next;
#	}
#	else {
#		$i++;
#		next;
#	}
}
$i=0;

#######################################################
# create subroutine to parse FASTA header
sub parse_header {

my $line = $_[0];
my @line;
my $orfline;
my @orfline;
my $start;
my $stop;
my $length;
my $orfotu;
my @orfotu;
my $orf;
my $otu;

		@line = split(/ /, $line); #process header line
			# 0 - >lcl|Otu1:2-310
			# 1 - ORF1_Otu1:1:309
		$orfline = $line[1];

		@orfline = split(/:/, $orfline);
			# 0 - ORF1_Otu1
			# 1 - 1
			# 2 - 309
		$start = $orfline[1];
		$stop = $orfline[2];
		$length = $stop - $start + 1;
		$orfotu = $orfline[0];

		@orfotu = split(/_/, $orfotu);
			# 0 - ORF1
			# 1 - Otu1
		$orf = $orfotu[0];
		$otu = $orfotu[1];

		return ($orf, $otu, $length);
 }


#!/home/terri/miniconda3/envs/myenv/bin/perl
# Dec. 13, 2019 by Teresita M. Porter
# Edit to handle otu id's formatted like: amplicon_otuID
# Script to grab longest CDS from ORFfinder nucleotide output (outfmt 1)
## NOT using cutoff = minimum number of base pairs required to retain CDS ##
# print to STDOUT here, in snakemake, redirect to an outfile
# USAGE perl get_longest_orf.plx orf_nt.fasta > cds.fasta.tmp

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
		$seq = $seq{$otu}{$longest};
		print STDOUT ">$otu\n$seq\n";
		$i++;
		next;
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
my $amplicon;
my $otu;
my $ampliconOtu;

		@line = split(/ /, $line); #process header line
			# 0 - >lcl|Amplicon_Otu1:2-310
			# 1 - ORF1_Amplicon_Otu1:1:309
		$orfline = $line[1];

		@orfline = split(/:/, $orfline);
			# 0 - ORF1_Amplicon_Otu1
			# 1 - 1
			# 2 - 309
		$start = $orfline[1];
		$stop = $orfline[2];
		$length = $stop - $start + 1;
		$orfotu = $orfline[0];

		@orfotu = split(/_/, $orfotu);
			# 0 - ORF1
			# 1 - Amplicon
			# 2 - Otu1
		$orf = $orfotu[0];
		$amplicon = $orfotu[1];
		$otu = $orfotu[2];
		$ampliconOtu = $amplicon."_".$otu;

		return ($orf, $ampliconOtu, $length);
 }


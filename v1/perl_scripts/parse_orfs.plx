#!/home/terri/miniconda3/envs/myenv/bin/perl
# Dec. 13, 2019 by Teresita M. Porter
# Script to grab matching orfs from nt and aa files
# calculate cutoffs to remove nt orfs with length outliers
# output matching filtered nt and aa orfs
# USAGE perl parse_orfs.plx orf_nt.fasta orf_aa.fasta

use strict;
use warnings;
use Data::Dumper;

# declare var
my $otu;
my $orf;
my $length;
my $ntseq;
my $aaseq;
my $lq = 25; #lower quartile
my $uq = 75; #upper quartile
my $percentile25;
my $percentile75;
my $iqr;
my $longest;
my $min;
my $max;
my $i=0;
my $outfile1 = $ARGV[0].".filtered";
my $outfile2 = $ARGV[1].".filtered";


# declare array
my @nt;
my @aa;
my @lengths;
my @otus;
my @longest;

# declare hash
my %ntLength; # key1 = otu, key2 = orf, value = length
my %ntSeq;
my %aaLength;
my %aaSeq;

my %match; # key1 = otu, key2 = orf, value = length

# read in files from ORFfinder
open (IN, "<", $ARGV[0]) || die "Can't open orf_nt.fasta: $!\n";
@nt = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Can't open orf_aa.fasta: $!\n";
@aa = <IN2>;
close IN2;

# hash length and seq from orf_nt.fasta
my ($ref1, $ref2) = parse_nt_orf(\@nt,\%ntLength, \%ntSeq);
%ntLength = %{$ref1};
%ntSeq = %{$ref2};

# hash length and seq from orf_aa.fasta
($ref1, $ref2) = parse_aa_orf(\@aa, \%aaLength, \%aaSeq);
%aaLength = %{$ref1};
%aaSeq = %{$ref2};

# keep records that match between nt and aa hashes
foreach $otu ( keys %ntLength) {
	foreach $orf (keys %{ $ntLength{$otu} }) {
		if (exists $aaLength{$otu}{$orf}) {
			$length = $ntLength{$otu}{$orf};
			$match{$otu}{$orf} = $length;
			push(@lengths, $length);
		}
	}
}

# calculate percentiles
$percentile25 = get_percentile(\@lengths, $lq);
$percentile75 = get_percentile(\@lengths, $uq);

# calculate interquartile range (IQR)
$iqr = $percentile75 - $percentile25;

# calculate lower cutoff: 25th percentile - 1.5*IQR
$min = $percentile25 - ($iqr * 1.5);

# calculate upper cutoff: 75th percentile + 1.5*IQR
$max = $percentile75 + ($iqr * 1.5);

# create outfiles
open (OUT1, ">>", $outfile1) || die "Cannot open orf_nt.fasta.filtered:$!\n";
open (OUT2, ">>", $outfile2) || die "Cannot open orf_aa.fasta.filtered:$!\n"; 

foreach $otu ( keys %match ) {
#	foreach $orf ( keys %{ $match{$otu} } ) {
		# for each otu, sorth the orfs by ascending length
		@longest =  sort { $match{$otu}{$a} <=> $match{$otu}{$b} } keys %{ $match{$otu} };
		$longest = $longest[-1]; #inner key to longest is at bottom of array
		$length = $match{$otu}{$longest};

		if ($length >= $min && $length <= $max) {
			$ntseq = $ntSeq{$otu}{$longest};	
			print OUT1 ">".$otu."\n".$ntseq."\n";
			$aaseq = $aaSeq{$otu}{$longest};
			print OUT2 ">".$otu."\n".$aaseq."\n";
		}
#	}
}

#######################################################
# create subroutine to parse FASTA header
sub parse_nt_header {

my $line = $_[0];
my @line;
my $orfline;
my @orfline;
my $start;
my $stop;
my $orfotu;
my @orfotu;
my $amplicon;
my $ampliconOtu;

		@line = split(/ /, $line); #process headier line
			# 0 - >lcl|Amplicon_Otu1:2:310
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

#######################################################
# parse out lengths and seqs from orf_nt.fasta

sub parse_nt_orf {

my @nt = @{$_[0]}; # access array reference
my %ntLength = %{$_[1]};
my %ntSeq = %{$_[2]};
my $i = 0;
my $line;
my $flag = 0;
my $orf;
my $otu;
my $length;
my $ntseq;

	while ($nt[$i]) {
		$line = $nt[$i];
		chomp $line;

		if ($flag==0 && $line =~/^>/) {
			($orf, $otu, $length) = parse_nt_header($line);
			$ntLength{$otu}{$orf} = $length;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag==1) {
			$ntseq = $line;
			$flag = 2;
			$i++;
			next;
		}
		elsif ($flag==2 && $line !~ /^>/) {
			$ntseq = $ntseq.$line;
			$i++;
			next;
		}
		else {
			$flag = 0;
			$ntSeq{$otu}{$orf} = $ntseq;
			next;
		}
	}
	$i=0;

	# add last seq to hash
	$ntSeq{$otu}{$orf} = $ntseq;

	return(\%ntLength, \%ntSeq);

}

#######################################################
# create subroutine to parse FASTA header
sub parse_aa_header {

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
		# 0 - >lcl|ORF1_Amplicon_Otu1:2:310
		# 1 - unnamed
		# 2 - protein
		# 3 - product,
		#4 - partial
		$orfline = $line[0];

		@orfline = split(/:/, $orfline);
			# 0 - >lcl|ORF1_Amplicon_Otu1
			# 1 - 2
			# 2 - 310
		$start = $orfline[1];
		$stop = $orfline[2];
		$length = $stop - $start + 1;
		$orfotu = $orfline[0];
		$orfotu =~ s/^>lcl\|//g;

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

#######################################################
# parse out lengths and seqs from orf_aa.fasta

sub parse_aa_orf {

my @aa = @{$_[0]}; # access array reference
my %aaLength = %{$_[1]};
my %aaSeq = %{$_[2]};
my $i = 0;
my $line;
my $flag = 0;
my $orf;
my $otu;
my $length;
my $aaseq;

	while ($aa[$i]) {
		$line = $aa[$i];
		chomp $line;

		if ($flag==0 && $line =~/^>/) {
			($orf, $otu, $length) = parse_aa_header($line);
			$aaLength{$otu}{$orf} = $length;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag==1) {
			$aaseq = $line;
			$flag = 2;
			$i++;
			next;
		}
		elsif ($flag==2 && $line !~ /^>/) {
			$aaseq = $aaseq.$line;
			$i++;
			next;
		}
		else {
			$flag = 0;
			$aaSeq{$otu}{$orf} = $aaseq;
			next;
		}
	}
	$i=0;

	# add last seq to hash
	$aaSeq{$otu}{$orf} = $aaseq;

	return(\%aaLength, \%aaSeq);

}

#######################################################
# Calculate percentile

sub get_percentile {

my @lengths = @{$_[0]};
my $percentile = $_[1];

my @sorted_lengths;
my $array_length;
my $decimal = $percentile/100;
my $index;
my $value;

	@sorted_lengths = sort { $a <=> $b } @lengths;
	$array_length = scalar(@sorted_lengths);
	$index = int($array_length * $decimal);
	$value = $sorted_lengths[$index];
	return($value);

}

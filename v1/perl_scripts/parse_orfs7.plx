#!/home/terri/miniconda3/envs/myenv/bin/perl
# Feb. 6,2020 by Teresita M. Porter
# Already translated, just keep longest, no need to match with a nt file
# calculate cutoffs to remove nt orfs with length outliers
# USAGE perl parse_orfs7.plx GenBank_pseudogene.nt.fasta 

use strict;
use warnings;
use Data::Dumper;

# declare var
my $accession;
my $orf;
my $length;
my $ntseq;
my $longest;
my $min;
my $max;
my $i=0;
my $outfile1 = $ARGV[0].".filtered";

# declare array
my @nt;
my @lengths;
my @otus;
my @longest;

# declare hash
my %ntLength; # key1 = otu, key2 = orf, value = length
my %ntSeq;

open (IN2, "<", $ARGV[0]) || die "Can't open orf_nt.fasta: $!\n";
@nt = <IN2>;
close IN2;

# hash length and seq from orf_nt.fasta
my ($ref1, $ref2) = parse_nt_orf(\@nt, \%ntLength, \%ntSeq);
%ntLength = %{$ref1};
%ntSeq = %{$ref2};

# print to file
open (OUT1, ">>", $outfile1) || die "Cannot open orf_nt.fasta.filtered:$!\n";

foreach $accession ( keys %ntLength ) {
	@longest =  sort { $ntLength{$accession}{$a} <=> $ntLength{$accession}{$b} } keys %{ $ntLength{$accession} };
	$longest = $longest[-1]; #inner key to longest is at bottom of array
	$ntseq = $ntSeq{$accession}{$longest};
	print OUT1 ">".$accession."\n".$ntseq."\n";
}

close OUT1;

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
my $accession;
my $length;
my $ntseq="";

	while ($nt[$i]) {
		$line = $nt[$i];
		chomp $line;

		if ($flag==0 && $line =~/^>/) {
			$ntseq="";
			($orf, $accession, $length) = parse_nt_header($line);
			$ntLength{$accession}{$orf} = $length;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag==1 && $line !~ /^>/) {
			$ntseq = $ntseq.$line;
			$i++;
			next;
		}
		elsif ($flag==1 && $line =~ /^>/) {
			$ntSeq{$accession}{$orf} = $ntseq;
			$flag=0;
			next;
		}
	}
	$i=0;

	# add last seq to hash
	$ntSeq{$accession}{$orf} = $ntseq;

	# testing 
	my $scalarSeq = keys (%ntSeq);
	my $scalarLength = keys (%ntLength);

	return(\%ntLength, \%ntSeq);

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
my $length;
my $orfacc;
my @orfacc;
my $orf;
my $accession;

		@line = split(/ /, $line); #process header line
		# 0 - >lcl|ABZ7793:1-117
		# 1 - ORF1_ABZ7793:0:116
		$orfline = $line[1];

		@orfline = split(/:/, $orfline);
			# 0 - ORF1_ABZ7793
			# 1 - 0
			# 2 - 116
		$start = $orfline[1];
		$stop = $orfline[2];
		$length = $stop - $start + 1;
		$orfacc = $orfline[0];
#		$orfacc =~ s/^>\w{3}\|//g;

		@orfacc = split(/_/, $orfacc);
			# 0 - ORF1
			# 1 - ABZ7793
		$orf = $orfacc[0];
		$accession = $orfacc[1];

		return ($orf, $accession, $length);
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

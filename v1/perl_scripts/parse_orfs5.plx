#!/home/terri/miniconda3/envs/myenv/bin/perl
# Feb. 5,2020 by Teresita M. Porter
# Already translated, just keep longest, no need to match with a nt file
# calculate cutoffs to remove nt orfs with length outliers
# output matching filtered nt and aa orfs
# USAGE perl parse_orfs5.plx GenBank_pseudogene.aa.fasta 

use strict;
use warnings;
use Data::Dumper;

# declare var
my $accession;
my $orf;
my $length;
my $aaseq;
my $longest;
my $min;
my $max;
my $i=0;
my $outfile1 = $ARGV[0].".filtered";

# declare array
my @aa;
my @lengths;
my @otus;
my @longest;

# declare hash
my %aaLength; # key1 = otu, key2 = orf, value = length
my %aaSeq;

open (IN2, "<", $ARGV[0]) || die "Can't open orf_aa.fasta: $!\n";
@aa = <IN2>;
close IN2;

# hash length and seq from orf_aa.fasta
my ($ref1, $ref2) = parse_aa_orf(\@aa, \%aaLength, \%aaSeq);
%aaLength = %{$ref1};
#print Dumper(\%aaLength);
%aaSeq = %{$ref2};

# print to file
open (OUT1, ">>", $outfile1) || die "Cannot open orf_aa.fasta.filtered:$!\n";

foreach $accession ( keys %aaLength ) {
	@longest =  sort { $aaLength{$accession}{$a} <=> $aaLength{$accession}{$b} } keys %{ $aaLength{$accession} };
	$longest = $longest[-1]; #inner key to longest is at bottom of array
	$aaseq = $aaSeq{$accession}{$longest};
	print OUT1 ">".$accession."\n".$aaseq."\n";
}

close OUT1;

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
my $orfacc;
my @orfacc;
my $orf;
my $accession;

		@line = split(/ /, $line); #process header line
		# 0 - >lcl|ORF6_AF301593.1:387:917
		# 1 - unnamed
		# 2 - protein
		# 3 - product,
		#4 - partial
		$orfline = $line[0];

		@orfline = split(/:/, $orfline);
			# 0 - >lcl|ORF6_AF301593.1
			# 1 - 387
			# 2 - 917
		$start = $orfline[1];
		$stop = $orfline[2];
		$length = $stop - $start + 1;
		$orfacc = $orfline[0];
		$orfacc =~ s/^>\w{3}\|//g;

		@orfacc = split(/_/, $orfacc);
			# 0 - ORF6
			# 1 - AF301593.1
		$orf = $orfacc[0];
#		print $orf."\n"; # testing
		$accession = $orfacc[1];

		return ($orf, $accession, $length);
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
my $accession;
my $length;
my $aaseq="";

	while ($aa[$i]) {
		$line = $aa[$i];
		chomp $line;

		if ($flag==0 && $line =~/^>/) {
			$aaseq="";
			($orf, $accession, $length) = parse_aa_header($line);
			$aaLength{$accession}{$orf} = $length;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag==1 && $line !~ /^>/) {
			$aaseq = $aaseq.$line;
			$i++;
			next;
		}
		elsif ($flag==1 && $line =~ /^>/) {
			$aaSeq{$accession}{$orf} = $aaseq;
			$flag=0;
			next;
		}
	}
	$i=0;

	# add last seq to hash
	$aaSeq{$accession}{$orf} = $aaseq;

	# testing 
	my $scalarSeq = keys (%aaSeq);
	my $scalarLength = keys (%aaLength);
#	print "scalarSeq $scalarSeq\t scalarLength $scalarLength\n";

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

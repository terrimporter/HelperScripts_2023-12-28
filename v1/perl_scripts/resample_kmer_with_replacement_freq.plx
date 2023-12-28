#!/usr/bin/perl
#April 19, 2012 edited to sample from a set of different kmer length files according to a set of frequencies
#NEW usage perl resample_kmer_with_replacement_freq.plx kmer_30.fasta kmer_60.fasta kmer_100.fasta kmer_200.fasta
#
#April 9, 2012 edited to sample WITH REPLACEMENT and grab a nubmer of kmers instead of a proportion
#March 30, 2012 by Terri Porter
#Script to randomly resample kmers from kmer.fasta without replacement
#usage perl resample_kmer.plx kmer.fasta

use strict;
use warnings;
use List::Util 'shuffle';

#declare var
my $sample = 1000000;#####reset this value log10 (1K to 1million)
my $prop_30 = 0.45;#####reset these as needed
my $prop_60 = 0.35;
my $prop_100 = 0.10;
my $prop_200 = 0.10;
my $sample_30;
my $sample_60;
my $sample_100;
my $sample_200;
my $k=0;

#declare array
my @in_30;
my @in_60;
my @in_100;
my @in_200;

#declare hash
my %kmer;
my %kmer_30;
my %kmer_60;
my %kmer_100;
my %kmer_200;
my %lines;

print "\nReading infiles...\n\n";
open (IN,"<",$ARGV[0]) || die "Error cannot read from kmer.fasta: $!\n";
@in_30 = <IN>;
close IN;
print "First of four files read.\n";

open (IN2,"<",$ARGV[1]) || die "Error cannot read from 2nd kmer.fasta: $!\n";
@in_60 = <IN2>;
close IN2;
print "Second of four files read.\n";

open (IN3,"<",$ARGV[2]) || die "Error cannot read from 3rd kmer.fasta: $!\n";
@in_100 = <IN3>;
close IN3;
print "Third of four files read.\n";

open (IN4,"<",$ARGV[3]) || die "Error cannot read from 4th kmer.fasta: $!\n";
@in_200 = <IN4>;
close IN4;
print "All four files read.\n";

#put all kmer seqs into hashes
print "\nCreating hash tables...\n\n";
print "Hashing for 30 bp kmer in progress...\n";
makeHash(\@in_30);
%kmer_30 = %kmer;
%kmer=();
print "Hash table for 30 bp kmer done.\n";

print "Hashing for 60 bp kmer in progress...\n";
makeHash(\@in_60);
%kmer_60 = %kmer;
%kmer=();
print "Hash table for 60 bp kmer done.\n";

print "Hashing for 100 bp kmer in progress...\n";
makeHash(\@in_100);
%kmer_100 = %kmer;
%kmer=();
print "Hash table for 100 bp kmer done.\n";

print "Hashing for 200 bp kmer in progress...\n";
makeHash(\@in_200);
%kmer_200 = %kmer;
%kmer=();
print "Hash table for 200 bp kmer done.\n";

open (OUT,">>","kmer_rand.fasta") || die "Error cannot write to outfile: $!\n";

print "\nResampling kmers according to specified proportions...\n\n";
$sample_30 = int($prop_30*$sample);
resample($sample_30,\%kmer_30);
print "3 0bp kmers resampling done.\n";

$sample_60 = int($prop_60*$sample);
resample($sample_60,\%kmer_60);
print "60 bp kmer resampling done.\n";

$sample_100 = int($prop_100*$sample);
resample($sample_100,\%kmer_100);
print "100 bp kmer resampling done.\n";

$sample_200 = int($prop_200*$sample);
resample($sample_200,\%kmer_200);
print "200 bp kmer resampling done.\n";

close OUT;

####################

sub makeHash {

#declare var
my $i=0;
my $line;
my $kmerID;
my $j;
my $line2;
my $kmerSeq;
my $kmerIDSeq;
#my $in = @_; #read array ref

#declare array
my @in = @{$_[0]};#dereference entire array

	while($in[$i]) {
		$line = $in[$i];
		chomp $line;

		if ($line =~ /^>/) {
			$line =~ s/^>//;
			$kmerID = $line;

			$j=$i+1;
			$line2 = $in[$j];
			chomp $line2;
			$kmerSeq = $line2;

			$kmerIDSeq = $kmerID."|".$kmerSeq;	
			$kmer{$kmerIDSeq} = 1;
		}
		$i++;
		$line=();
		$kmerID=();
		$j=();
		$line2=();
		$kmerSeq=();
		$kmerIDSeq=();
	}
	$i=0;
	return %kmer;

}

####################

sub resample {
#RESAMPLE WITH REPLACEMENT

#declare var
my $count;
my $range;
my $scalar=0;
my $number = shift;
my $random_index;
my $kmerIDSeq;
my $kmerID;
my $kmerSeq;

#declare array
my @keys;
my @lines;
my @kmerIDSeq;

#declare hash
my %kmer = %{$_[0]};

	@keys = keys(%kmer); #put hash keys into array so they can be sampled using the random index values 
	$count = keys(%kmer);
	$range = $count+1;

	while($scalar < $number) {
		$random_index = int(rand($range));
		push (@lines, $random_index);
		$scalar = scalar(@lines);
		$kmerIDSeq = $keys[$random_index];
		@kmerIDSeq = split(/\|/,$kmerIDSeq);
		$kmerID = $kmerIDSeq[0];
		$kmerSeq = $kmerIDSeq[1];
		print OUT ">$k|$kmerID\n$kmerSeq\n";
		
		$k++;
		$random_index=();
		$kmerIDSeq=();
		@kmerIDSeq=();
		$kmerID=();
		$kmerSeq=();
	}
	
	$scalar=0;
	@lines=();
	@keys=();
	$count=();
	$range=();

}

####################


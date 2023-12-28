#!/usr/bin/perl
#June 15, 2011 by Terri Porter
#Script to create sampling curve from velvet data, afg.parsed
#usage perl resample_reads.plx afg.parsed

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $contig;
my $aref; 
my $key;
my $element;
my $num_unique_reads;
my $num_reads;
my $num_contigs;
my $contents;
my $j=1;
my $random_read_index;
my $random_read;
my $num_intervals;
my $fraction;
my $current_fraction;
my $rounded;
my $k=1;
my $num_to_sample;
my $random_reads;
my $unique_reads;
my $l=0;
my $retrieved_contig;
my $read;
my $readid;

#declare array
my @afg;
my @line;
my @aref;
my @all_reads;
my @unique_reads;
my @random_reads;

#declare hash
#my %contig_contents;
my %all_reads;
my %random_reads;
my %retrieved_contigs;
my %read_contigs;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@afg = <IN>;
close IN;

open (OUT,">>","table.txt") || die ("Error cannot write to outfile: $!\n");

#create hash with keys as contigIDs and values as arrays of readIDs
#while ($afg[$i]) {
#	$line = $afg[$i];
#	chomp $line;

#	@line = split(/ /,$line);
#	$contig = shift (@line);
#	$contig_contents{$contig} = [@line];	

	#create list of all reads
#	foreach (@line) {
#		push (@all_reads, $_)
#	}

#	$i++;

#}
#print "Got contig read hash\n";

#instead of hash with keys as contigIDs and values as arrays of readIDs USE
#hash with keys as readIDs and values as arrays of contigs to make searching MUCH faster

while ($afg[$i]) {
	$line = $afg[$i];
	chomp $line;

	@line = split(/ /,$line);
	$contig = shift(@line);

	foreach (@line) {
		$readid = $_;
		push @{$read_contigs{$readid}},$contig;
		push (@all_reads, $readid);
	}
	$i++;
}

#get list of unique reads only
%all_reads = map { $_, 1} @all_reads;
@unique_reads = keys %all_reads;

$num_unique_reads = scalar(@unique_reads);
print "Unique reads: $num_unique_reads\n"; #test

get_interval();



####################get interval and number of reads to sample

sub get_interval {

#print "Please enter number of intervals (ex. 10, 100, 1000):\n";
#$num_intervals = <STDIN>;
#chomp $num_intervals;
$num_intervals=10;#hard code so can run w/o user intervention with qsub
$fraction = $num_unique_reads/$num_intervals;

while ($k <= ($num_intervals)) {

	$current_fraction = $fraction*$k;
	$num_to_sample = int($current_fraction + 0.5);	
	print "Num_to_sample: $num_to_sample\n";#test
	print "Picking randomly from array...";
	pick_randomly_from_array();
	print "random reads chosen for interval $k\n";
	$num_reads = keys(%random_reads);#test
	print "reads: $num_reads\t";#test
	print "finding contigs...";
	count_contigs();
	$num_contigs = keys %retrieved_contigs;
	print "contigs: $num_contigs\n";
	print OUT "$k\t$num_reads\t$num_contigs\n";
	for (keys %random_reads) {
		delete $random_reads{$_};
	}
	for (keys %retrieved_contigs) {
		delete $retrieved_contigs{$_};
	}
	$k++;
}
close OUT;
}

####################figure out contigs, and count them

sub count_contigs {

while (my ($key,$value) = each (%random_reads)) {
	$random_read = $key;
#	print "for each random key..."; 
	
	while (my ($key2,$aref) = each %read_contigs) {
		if ($random_read eq $key2) {
			@aref = @$aref;	
			foreach (@aref) {
				$retrieved_contig = $_;
				$retrieved_contigs{$retrieved_contig} = 1;
			}
		}
	}
			
#	while (my ($key2,$aref) = each %contig_contents) {
#		$retrieved_contig = $key2;
#		@aref = @$aref;
#		print "checking for contigs\n";
#		while ($aref[$l]) {
#			$read = $aref[$l];
#			if ($read =~ /$random_read/g) {
#				$retrieved_contigs{$retrieved_contig} = 1;
#			}
#			$l++;
#		}
#		$l=0;
#		$retrieved_contig=();
#		@aref=();
#	}
}
print "done.\n";
}

####################pick random elements from array without repeats/replacement

sub pick_randomly_from_array {

while ($j <= $num_to_sample) {
	$random_read = $unique_reads[rand @unique_reads];
#	print "$random_read_index\n"; #test
#	$random_read = $unique_reads[$random_read_index];
#	print "random_read: $random_read\n";#test
	redo if defined $random_reads{$random_read};
	$random_reads{$random_read} = 1;
	$j++;
}
$j=1;
print "done.\n";
}

####################check that hash looks ok
sub check_hash {

while (my($key,$aref) = each %read_contigs) {	
	@aref = @$aref;
	print "$key => @aref\n";
}

}
####################

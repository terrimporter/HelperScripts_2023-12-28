#!/usr/bin/perl
#Nov. 8, 2013 by Terri Porter
#Written for Shadi Shokralla, Guelph
#Script to create all combinations of k-mer oligo
#From this list, only print the ones with an edit distance > x

use strict;
use warnings;
use Algorithm::Combinatorics qw(variations_with_repetition);
use Text::Levenshtein qw(distance);

#declare var
my $length;
my $distance;
my $iter;
my $combo;
my $i=0;
my $scalar1=1;
my $scalar2=0;
my $x;
my $oligo;
my $first;
my $flag;

#declare array
my @bases= ("A","C","G","T");
my @oligos;
my @distances;
my @temp;
my @temp2;
my @sorted;

print STDOUT "Enter oligo length (bp):\n";
$length = <STDIN>;
chomp $length;

print STDOUT "Press 1 = HAMMING DISTANCE or 2 = EDIT DISTANCE:\n";
$flag = <STDIN>;
chomp $flag;

print STDOUT "Enter distance cutoff:\n";
$distance = <STDIN>;
chomp $distance;

$iter = variations_with_repetition(\@bases, $length);

while ($combo = $iter->next) {
	$oligo = join('',@{ $combo });
	push(@oligos,$oligo);
}

while ($scalar1 > $scalar2) {
	$scalar1 = scalar(@oligos);
	$first = shift @oligos; 
	@temp = @oligos; 

	while ($temp[$i]) {
		$oligo = $temp[$i];

		if ($flag == 1) {
			$x = hamming($first,$oligo);
		}

		elsif ($flag == 2) {
			$x = editdistance($first,$oligo);
		}

		if ($x > $distance) {
			push (@temp2,$oligo);
		}
		$i++;
		$oligo=();
		$x=();
	}
	$i=0;
	@temp=();
	push (@temp2,$first);
	@oligos=();
	@oligos = @temp2;
	@temp2=();
	$scalar2 = scalar(@oligos);
}

if ($flag == 1 ) {
	print "Final oligos with HAMMING DISTANCE > $distance:\n";
}
elsif ($flag == 2) {
	print "Final oligos with EDIT DISTANCE > $distance:\n";
}

@sorted = sort(@oligos);
foreach $oligo (@sorted) {
	print "$oligo\n";
}

####################

#copied revised subroutine from http://www.perlmonks.org/?node_id=500235

sub hamming {

my $mismatch=0;
my $j=0;

return 0 if $first eq $oligo;

for ($j=0; $j<$length; $j++) {
	++$mismatch if substr($first,$j,1) ne substr($oligo,$j,1);
}
return $mismatch;

}

####################

sub editdistance {

my $editdistance = 0;

$editdistance = distance($first,$oligo);
return $editdistance;

}

####################

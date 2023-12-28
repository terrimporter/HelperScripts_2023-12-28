#!/usr/bin/perl
# Oct. 28/21 by Teresita M. Porter
# Script to combine the results from RDP Classifier LOOT partial files, focusing on the first table
# USAGE perl combine_loot_parts.plx

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils 'pairwise';

# declare vars
my $dir = "/home/tep001/lv_workspace/UNITE_ITSClassifier/v2/mydata/testing";
my $i=0;
my $filename;
my $phylum;
my $species;
my $genus;
my $root;
my $family;
my $class;
my $kingdom;
my $order;
my $val;
my $rank;
my $last;
my $outfile = "summed_table.csv";

# declare arrays
my @files;
my @in;
my @vals;
my @vals_prev;
my @root;
my @kingdom;
my @phylum;
my @class;
my @order;
my @family;
my @genus;
my @species;

# declare hashes
my %phylum; #key = filename, value = string of table values, tab delimited
my %species;
my %genus;
my %root;
my %family;
my %class;
my %kingdom;
my %order;

# read in a directory of filenames
opendir(DIR, $dir) || die "Cannot open dir: $!\n";
@files = readdir(DIR);
closedir DIR;

# skip over files that start with a dot
# open file, grab table, hash contents by rank
while($files[$i]) {
	$filename = $files[$i];
	chomp $filename;

	# skip over filenames that start with a dot
	if ($filename =~ /^\./) {
		$i++;
		next;
	}
	else {
		open(IN, "<", $dir."/".$filename) || die "Cannot open filename: $!\n";
		@in = <IN>;
		close IN;

		# only need lines 12-19 in each file to grab the first table
		$phylum = $in[11];
		chomp $phylum;
		$phylum{$filename} = $phylum;

		$species = $in[12];
		chomp $species;
		$species{$filename} = $species;

		$genus = $in[13];
		chomp $genus;
		$genus{$filename} = $genus;
	
		$root = $in[14];
		chomp $root;
		$root{$filename} = $root;

		$family = $in[15];
		chomp $family;
		$family{$filename} = $family;

		$class = $in[16];
		chomp $class;
		$class{$filename} = $class;

		$kingdom = $in[17];
		chomp $kingdom;
		$kingdom{$filename} = $kingdom;

		$order = $in[18];
		chomp $order;
		$order{$filename} = $order;

	}

	$i++;
}
$i=0;

# sum vals at each rank
$rank = 'root';
@root = sum_parts(\%root, $rank);

$rank = 'kingdom';
@kingdom = sum_parts(\%kingdom, $rank);

$rank = 'phylum';
@phylum = sum_parts(\%phylum, $rank);

$rank = 'class';
@class = sum_parts(\%class, $rank);

$rank = 'order';
@order = sum_parts(\%order, $rank);

$rank = 'family';
@family = sum_parts(\%family, $rank);

$rank = 'genus';
@genus = sum_parts(\%genus, $rank);

$rank = 'species';
@species = sum_parts(\%species, $rank);

# now put it all together in a proper table
open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

print OUT "Rank,100-95,,94-90,,89-80,,79-70,,69-60,,59-50,,49-40,,39-30,,29-20,,19-10,,9-0,,CorrectAssts,TotalSeqs\n";
$root = join(",","Root",@root);
print OUT $root."\n";
$kingdom = join(",","Kingdom",@kingdom);
print OUT $kingdom."\n";
$phylum = join(",","Phylum",@phylum);
print OUT $phylum."\n";
$class = join(",","Class",@class);
print OUT $class."\n";
$order = join(",","Order",@order);
print OUT $order."\n";
$family = join(",","Family",@family);
print OUT $family."\n";
$genus = join(",","Genus",@genus);
print OUT $genus."\n";
$species = join(",","Species",@species);
print OUT $species."\n";

close OUT;

####################

sub sum_parts {

	# for each key in hash, pop off last val, sum current with previous to get an overall total
	@vals_prev = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
				  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	my %hash = %{$_[0]};
	my $rank = $_[1];

	foreach $val (values %hash) {
		@vals = split(/\t/,$val);
		# verify rank, and shift it off top of array
		$rank = shift(@vals);
		if ($rank eq $rank) {
			# pop last value off end of array
			$last = pop(@vals);
		
			# sum current with previous
			@vals_prev = pairwise { $a + $b } @vals_prev, @vals;
		}
		else {
			print "Expected root rank not found, instead found $rank\n";
		}
	}

	return @vals_prev;	
}

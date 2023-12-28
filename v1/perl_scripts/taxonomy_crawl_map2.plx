#!/usr/bin/perl
#April 17, 2012 by Terri Porter
#Script to map taxid.parsed (from taxonomy_crawl.plx) to readids using readid_taxonid_out-ex.txt (from MEGAN)
#usage perl taxonomy_crawl_map.plx readid_taxonid_out-ex.txt taxid.parsed

use warnings;
use strict;

#declare var
my $i=0;
my $line;
my $readid;
my $taxonid;
my $string;
my $new_string;
my $j=0;
my $count=0;
my $filename;

#declare arrays
my @map;
my @lineage;
my @line;
my @string;

#declare hash
my %map;
my %taxonidsSeen;

open (MAP,"<",$ARGV[0]) || die "Error cannot read mapping file: $!\n";
#@map = <MAP>;
#print "Finished reading merged-ex.txt\n";
#close MAP;

open (LINEAGE,"<",$ARGV[1]) || die "Error cannot read lineage file: $!\n";
@lineage = <LINEAGE>;
print "Finished reading taxid.parsed.sort.uniq\n";
close LINEAGE;

#parse mapping file into a hash indexed by taxonid with pipe-delimited readids
#while($map[$i]) {
print "Hashing...\n";
while(<MAP>) {
	$line = $_;
	chomp $line;
#	print "$line\n";

	@line = split(", ", $line);
	#print "array line: @line\n"; #test
	$readid = $line[0];
	$taxonid = $line[1];
#	print "$readid\t$taxonid\n";
	
	if ( exists $map{$taxonid} ) {
		$string = $map{$taxonid};
		$new_string = $string."|".$readid;
		$map{$taxonid} = $new_string;
#		print "$taxonid\t$new_string\n";
	}
	else {
		$map{$taxonid} = $readid;
#		print ".";
	}
	
#	$i++;
#	print ".";
	$line=();
	@line=();
	$readid=();
	$taxonid=();
	$string=();
	$new_string=();
}
#$i=0;
close MAP;
print "Finished hashing...\n";

#test hash
print "Counting hash keys\n";
while(my($key,$value) = each(%map)) {
#	print "$key => $value\n";
	$count++;
}
print "There were $count elements in map hash\n";

$filename = $ARGV[0].".parsed";

open (OUT,">>",$filename) || die "Error cannot write outfile: $!\n";

#for each taxonid\tlineage, prefix with readid
print "Mapping taxonid to lineage\n";
while ($lineage[$i]) {
	$line = $lineage[$i];
	chomp $line;
	
	@line = split(/\t/,$line);
#	print "array line: @line\n";
	$taxonid = $line[0];
#	$taxonidsSeen{$taxonid} = 1;
#	print "taxonid: $taxonid\n";

	if (exists $taxonidsSeen{$taxonid} ) {
		$i++;
		next;
	}
	else {
		$string = $map{$taxonid};
		@string = split(/\|/, $string);
#	print "array string: @string\n";

		while ($string[$j]) {
			$readid = $string[$j];
		
			print OUT "$readid\t$line\n";

			$j++;
			$readid=();
		}
		$taxonidsSeen{$taxonid} = 1;
		$j=0;
	}

	$i++;
	$line=();
	@line=();
	$taxonid=();
	$string=();
	@string=();
}
$i=0;

close OUT;

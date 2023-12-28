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

#declare arrays
my @map;
my @lineage;
my @line;
my @string;

#declare hash
my %map;
my %taxonidsSeen;

open (MAP,"<",$ARGV[0]) || die "Error cannot read mapping file: $!\n";
@map = <MAP>;
close MAP;

open (LINEAGE,"<",$ARGV[1]) || die "Error cannot read lineage file: $!\n";
@lineage = <LINEAGE>;
close LINEAGE;

#parse mapping file into a hash indexed by taxonid with pipe-delimited readids
while($map[$i]) {
	$line = $map[$i];
	chomp $line;
	
	@line = split(", ", $line);
	#print "array line: @line\n"; #test
	$readid = $line[0];
	$taxonid = $line[1];
	
	if ( exists $map{$taxonid} ) {
		$string = $map{$taxonid};
		$new_string = $string."|".$readid;
		$map{$taxonid} = $new_string;
	}
	else {
		$map{$taxonid} = $readid;
	}
	
	$i++;
	$line=();
	@line=();
	$readid=();
	$taxonid=();
	$string=();
	$new_string=();
}
$i=0;

#test hash
while(my($key,$value) = each(%map)) {
#	print "$key => $value\n";
	$count++;
}
print "There were $count elements in map hash\n";

open (OUT,">>","readid_taxonid.parsed") || die "Error cannot write outfile: $!\n";

#for each taxonid\tlineage, prefix with readid
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

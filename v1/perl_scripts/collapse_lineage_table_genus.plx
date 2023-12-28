#!/usr/bin/perl
#Oct.12, 2012 by Terri Porter
#Script to collapse entries in readid_taxonid.parsed.awk, where readid_taxonid.parsed had it's columsn reformatted by awk, where readid_taxonid.parsed came from taxonomy_crawl_map2.plx
#usage perl collapse_lineage_table.plx readid_taxonid.parsed.awk

use strict;
use warnings;

#declare var
my $line;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $readid;
my $taxonid;
my $readnum="";
my $lineage="";
my $i=0;
my $totalFields;
my $readidField;
my $taxonidField;
#declare array
my @line;
my @unsorted;
my @sorted;

#declare hash
my %lineage;

open (IN, "<", $ARGV[0]) || die "Error cannot read infile: $!\n";

#process the file line by line because it's huge!
print "Populating hash\n";
while(<IN>) {
	$line = $_;
	chomp $line;

	@line = split(' ',$line);
	$kingdom = $line[0];
	$phylum = $line[1];
	$class = $line[2];
	$order = $line[3];
	$family = $line[4];
	$genus = $line[5];
	$totalFields = scalar(@line);
	$readidField = $totalFields-1;
	$readid = $line[$readidField];
	$taxonidField = $totalFields-2;
	$taxonid = $line[$taxonidField];
#	print "$readnum\n";

	$lineage = $kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus;

	#if lineage already in hash, just update the readnum
	if (exists $lineage{$lineage}) {
		$readnum = $lineage{$lineage};
		$readnum += 1;
		$lineage{$lineage} = $readnum;
	}

	#otherwise, just add it to the hash
	else {
		$lineage{$lineage} = 1;
	}

	$line=();
	@line=();
	$kingdom=(); #column [0]
	$phylum=(); # [1]
	$class=(); # [2]
	$order=(); # [3]
	$family=(); #[4]
	$genus=(); # [5]
	$readid=();
	$taxonid=();
	$readnum="";
	$lineage="";

}
close IN;

#put hash contents into unsorted array
print "Sorting\n";
while (($lineage,$readnum) = each %lineage) {
	
	#to get a tab-separated entry line
	$line = $lineage."\t".$readnum; 

	#feed entries into unsorted array
	push(@unsorted,$line);

	$line=();
}

#Schwartzian transform
@sorted = map {$_->[0]}
	      sort { $a->[0] cmp $b->[0] || 
			     $a->[1] cmp $b->[1] || 
    			 $a->[2] cmp $b->[2] || 
	    		 $a->[3] cmp $b->[3] ||  
		    	 $a->[4] cmp $b->[4] ||
		 	     $a->[5] cmp $b->[5] }
	      map { [$_, split("\t")] }
	      @unsorted;

#print sorted data into the outfile

open (OUT, ">>", "genus.table.sorted") || die "Error cannot open outfile: $!\n";

while ($sorted[$i]) {
	$line = $sorted[$i];
	
	print OUT "$line\n";
	
	$i++;
	$line=();
}
$i=0;

close OUT;

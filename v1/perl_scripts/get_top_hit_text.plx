#!/usr/bin/perl
#Feb.8,2011 by Terri Porter
#Script to parse -outfmt 0 (pairwise text) and check top hit to see if accession matches query accession
#usage $perl get_top_hit_text.plx file.blast
#modify to get gi numbers for each hit, grab lineage info for these hits, compare these with the query

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

#var
my $line;
my $flag=0;
my $query;
my $i=0;
my $outfile;
my $outfile2;
my $hit;
my $factory;
my $file;
my $j=0;
my $gi;
my $gb;
my $factory2;
my $file2;
my $k=0;
my $l=0;
my $gi2;
my $genus;
my $species;
my $part;
my $outfile3;

#array
my @blast;
my @query;
my @line;
my @hit;
my @fasta_hits;
my @gi;
my @gb;
my @gb_file;
my @part;

#parse blast+ -outfmt 0 blast results to grab just top hit

open (IN,"<",$ARGV[0]) || die ("Error cannot read blast file: $!\n");
@blast = <IN>;
close IN;

$outfile = $ARGV[0];
chomp $outfile;
#$outfile2 = $outfile.".parsed";

#open (OUT,">>", $outfile2) || die ("Error cannot write to blast.parsed file: $!|n");

$outfile3 = $outfile.".parsed.all";
open (OUT2,">>",$outfile3) || die ("Error cannot write to blast.parsed.all file: $!\n");

while ($blast[$i]) {
	$line = $blast[$i];
	chomp $line;
	if ($flag==0) {
		if ($line =~ /^Query\=/) {
			$line =~ /Query\=\s+(\w+)\|\S+\|\S+\|\S+/g;
			$query = $1;
			push(@query,$query);
			$flag=1;
			print "Processing query: $query\n";#test
		}	
	}
	elsif ($flag==1) {
#		print "flag reset\n";#test
		if ($line =~ /^\w+\|\w+\.\d+\|/g) {
			#@line = split(/\|/,$line);
			#$hit=$line[1];
			#push(@hit,$hit);
#			print "Found hit line\n";#test
			if ($line =~ /(Uncultured|uncultured)/) {
				$i++;
				next;
			}
			elsif ($line =~ /sp\./) {
				$i++;
				next;
			}
			else  {
				@line = split(/\|/,$line);
				$hit=$line[1];
				push(@hit,$hit);
				print OUT2 "$query\t$line\n";
				$flag=0;
			}
			@line=();
		}

	}
	$i++;
}
close OUT;

#grab fasta files for each hit, parse gi number
my $test2 = scalar(@query);
my $test = scalar(@hit);
print "\nI have $test2 queries in query array\n";#test
print "\nI have $test hits in hit array\n";#test

$factory = Bio::DB::EUtilities -> new (	-eutil	=> 'efetch',
					-db	=> 'nucleotide',
					-id	=> \@hit,
					-rettype=> 'fasta');

$file = 'hits.fasta';
$factory -> get_Response(-file => $file);

open (IN,"<","hits.fasta") || die ("Error cannot read from hits.fasta: $!\n");
@fasta_hits = <IN>;
close IN;

open (OUT2,">>","gi.hits") || die ("Error cannot write to gi.hits: $!\n");

open (OUT3,">>","gb.hits") || die ("Error cannot write to gb.hits: $!\n");

while ($fasta_hits[$j]) {
	$line = $fasta_hits[$j];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>\w+\|(\d+)\|\w+\|(\S+)\|/;
		$gi = $1;
		push(@gi,$gi);
		print OUT2 "$gi\n";
		$gb = $2;
		push(@gb,$gb);
		print OUT3 "$gb\n";
	}
	$j++;
}
close OUT2;
close OUT3;

#use gi numbers to get genbank records

$factory2 = Bio::DB::EUtilities -> new (	-eutil	=> 'efetch',
						-db	=> 'nucleotide',
						-rettype=> 'genbank',
						-id	=> \@gi);

$file2 = 'hits.gb';
$factory2 -> get_Response (-file => $file2);

#parse gb records to get lineage (taxonomy)

open (IN,"<","hits.gb") || die ("Error cannot read hits.gb: $!\n");
@gb_file = <IN>;
close IN;

open (OUT,">>","hits.name") || die ("Error cannot write hits.name: $!\n");

while ($gb_file[$k]) {
	$line = $gb_file[$k];
	chomp $line;

	if ($line =~ /Seq-entry ::=/) {
		print "found seq-entry ::=\n";#test
		$a = 1;
		$k++;
		next;
	}
	elsif ($a==1) {
		if ($line =~ /genbank\s{1}\{/) {
			$a=2;
			$k++;
			next;
		}
	}
	elsif ($a==2) {
		if ($line =~ /gi\s{1}\d+\s{1}/) {
			$line =~ /gi\s{1}(\d+)\s{1}/;
			$gi2 = $1;
			print "found gi2\n";#test
			print OUT "\n$gi2\t";
			$a=3;
			$k++;
			next;
		}
	}
	elsif ($a==3) {
		if ($line =~ /title "/) {
			$line =~ /title "(.+)/g;
			$part = $1;
			@part = split(/ /,$part);
			$genus = $part[0];
			$species = $part[1];
			print OUT "$genus\t$species\t";
			$a=4;
			$k++;
			next;
		}
	}
	#elsif ($a==3) {
	#	if ($line =~ /binomial/) {
	#		$a=4;
	#		$k++;
	#		next;
	#	}
#	}
#	elsif ($a==4) {
#		if ($line =~ /genus "/) {
#			$line =~ /genus "(\w+)"/;
#			$genus = $1;
#			print "found genus\n";#test
#			print OUT "$genus\t";
#			$a=5;
#			$k++;
#			next;
#		}
#	}
#	elsif ($a==5) {
#		if ($line =~ /species "/) {
#			$line =~ /species "(.+)"/g;
#			$species = $1;
#			print "found species\n";#test
#			print OUT "\t$species";
#			$a=6;
#			$k++;
#			next;
#		}
#		else {
#			$a=6;
#			$k++;
#			next;
#		}
#	}
	$k++;
}
close OUT;

#use gi numbers to grab lineage info

#compare hit lineage with query lineage info

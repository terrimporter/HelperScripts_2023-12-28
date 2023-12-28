#!/usr/bin/perl
#Feb.8,2011 by Terri Porter
#Script to parse -outfmt 0 (pairwise text) and check top hit to see if accession matches query accession
#usage $perl get_top_hit_text.plx file.blast
#modify to get gi numbers for each hit, grab lineage info for these hits, compare these with the query
#modify to also compare names
#usage $perl NGS_pipeline_part3.plx file.blast name.query

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
my $outfile2b;
my $first;

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
print "Enter name of file.blast:\n";
$first = <STDIN>;
chomp $first;

open (IN,"<",$first) || die ("Error cannot read blast file: $!\n");
@blast = <IN>;
close IN;

$outfile = $first;
#$outfile2 = $outfile.".parsed";
$outfile2b = $outfile.".parsed.all";

#open (OUT,">>", $outfile2) || die ("Error cannot write to blast.parsed file: $!|n");
open (OUTb,">>",$outfile2b) || die ("Error: $!\n");

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
		if ($line =~ /^\w+\|\w+\.\d+\|/g) {
			if ($line =~ /(Uncultured|uncultured)/) {
				$i++;
				next;
			}
			elsif ($line =~ /sp\./) {
				$i++;
				next;
			}
			else {
				@line = split(/\|/,$line);
				$hit = $line[1];
				push(@hit,$hit);
				print OUTb "$query\t$line\n";
				$flag=0;
			}
			@line=();
		}
	}
	$i++;
}
#close OUT;
close OUTb;

#grab fasta files for each hit, parse gi number

#$factory = Bio::DB::EUtilities -> new (	-eutil	=> 'efetch',
#					-db	=> 'nucleotide',
#					-id	=> \@hit,
#					-rettype=> 'fasta');

#$file = 'hits.fasta';
#$factory -> get_Response(-file => $file);

#open (IN,"<","hits.fasta") || die ("Error cannot read from hits.fasta: $!\n");
#@fasta_hits = <IN>;
#close IN;

#open (OUT2,">>","gi.hits") || die ("Error cannot write to gi.hits: $!\n");

#open (OUT3,">>","gb.hits") || die ("Error cannot write to gb.hits: $!\n");

#while ($fasta_hits[$j]) {
#	$line = $fasta_hits[$j];
#	chomp $line;

#	if ($line =~ /^>/) {
#		$line =~ /^>\w+\|(\d+)\|\w+\|(\S+)\|/;
#		$gi = $1;
#		push(@gi,$gi);
#		print OUT2 "$gi\n";
#		$gb = $2;
#		push(@gb,$gb);
#		print OUT3 "$gb\n";
#	}
#	$j++;
#}
#close OUT2;
#close OUT3;

#use gi numbers to get genbank records

#$factory2 = Bio::DB::EUtilities -> new (	-eutil	=> 'efetch',
#						-db	=> 'nucleotide',
#						-rettype=> 'genbank',
#						-id	=> \@gi);

#$file2 = 'hits.gb';
#$factory2 -> get_Response (-file => $file2);

#parse gb records to get lineage (taxonomy)

#open (IN,"<","hits.gb") || die ("Error cannot read hits.gb: $!\n");
#@gb_file = <IN>;
#close IN;

#open (OUT,">>","hits.name") || die ("Error cannot write hits.name: $!\n");

#while ($gb_file[$k]) {
#	$line = $gb_file[$k];
#	chomp $line;

#	if ($line =~ /Seq-entry ::=/) {
#		print "found seq-entry ::=\n";#test
#		$a = 1;
#		$k++;
#		next;
#	}
#	elsif ($a==1) {
#		if ($line =~ /genbank\s{1}\{/) {
#			$a=2;
#			$k++;
#			next;
#		}
#	}
#	elsif ($a==2) {
#		if ($line =~ /gi\s{1}\d+\s{1}/) {
#			$line =~ /gi\s{1}(\d+)\s{1}/;
#			$gi2 = $1;
#			print "found gi2\n";#test
#			print OUT "\n$gi2\t";
#			$a=3;
#			$k++;
#			next;
#		}
#	}
#	elsif ($a==3) {
#		if ($line =~ /title "/) {
#			$line =~ /title "(.+)/g;
#			$part = $1;
#			@part = split(/ /,$part);
#			$genus = $part[0];
#			$species = $part[1];
#			print OUT "$genus\t$species\t";
#			$a=4;
#			$k++;
#			next;
#		}
#	}
#	$k++;
#}
#close OUT;

#compare names

#var
my $filename_query;
my $filename_parsed_all;
my $m=0;
my $query_acc;
my $hit_acc;
my $hit_genus;
my $hit_species;
my $exact_match_counter=0;
my $n=0;
my $line2;
my $query_genus;
my $query_species;
my $species_match_counter=0;
my $genus_match_counter=0;

#array
my @query2;
my @blast2;
my @line2;
my @output;

@output = qx(ls -lhrt);
foreach $line (@output) {
	print $line;
}

print "Enter path/filename to name.query (ex. /home/terri/ITS_NGS/name.query):\n";
$filename_query = <STDIN>;
chomp $filename_query;

open (QUERY2,"<",$filename_query) || die ("Error cannot open query name file: $!\n");
@query2=<QUERY2>;
close QUERY2;

open (BLAST2,"<",$outfile2b) || die ("Error cannot open hit name file: $!\n");
@blast2=<BLAST2>;
close BLAST2;

while ($blast2[$m]) {
	$line = $blast2[$m];
	chomp $line;
	if ($line =~ /^\w+\t\w+\|\w+\.\d+\|\s+\w+\s+\S+\s+/) {
		$line =~ /^(\w+)\t\w+\|(\w+)\.\d+\|\s+(\w+)\s+(\S+)\s+/;
		$query_acc = $1;
		$hit_acc = $2;
		$hit_genus = $3;
		$hit_species = $4;
		#print "found line pattern\n";#test

		#check for exact sequence match
		if ($query_acc eq $hit_acc) {
			$exact_match_counter++;
		}

		#check for species match and genus match
		while ($query2[$n]) {
			$line2 = $query2[$n];
			chomp $line2;
			if ($line2 =~ /^$query_acc/g) {
				@line2 = split(/\t/,$line2);
				$query_genus = $line2[1];
				$query_species = $line2[2];
				print "query: $query_species\thit: $hit_species\n";#test
				#if ($query_species eq $hit_species) {
					#print "query: $query_species\thit: $hit_species\n";#test
					#$species_match_counter++;
				#}
				if ($query_genus eq $hit_genus) {
					$genus_match_counter++;
					if ($query_species eq $hit_species) {
						#if ($query_species ne "sp.") {
							$species_match_counter++;
						#}
					}
				}
			}
			$n++;
		}
	}
	$n=0;
	$m++;
}
print "\nexact matches = $exact_match_counter\nspecies matches = $species_match_counter\ngenus matches = $genus_match_counter\n";


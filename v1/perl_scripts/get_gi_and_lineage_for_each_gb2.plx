#!/usr/bin/perl
#March 31, 2011 by Terri Porter
#Script to use a list of full gb numbers to get gi numbers
#usage $ perl get_gi_for_each_gb.plx hit_gb.list
#May 11, 2011 modified to parse query_gi|query|gb\thit_gb format

use strict;
use warnings;
use Bio::DB::EUtilities;

#declare var
my $i=0;
my $line;
my $factory;
my $one_record;
my $j=0;
my $hit_gb;
my $hit_gi;
my $flag=0;
my $lineage;
my $query_map;
my $k=0;
my $l=0;
my $species;
my $gb;

#declare array
my @in;
my @gb;
my @hit_gis;
my @organism;
my @read;
my @line;
my @query_map;
my @read2;
my @in3;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@in= <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$query_map = $line[0];
	push(@query_map,$query_map);
	$hit_gb = $line[1];
	push (@gb,$hit_gb);
	$i++;
}

$factory = Bio::DB::EUtilities-> new(	-eutil	=> 'efetch',
					-db	=>'nucleotide',
					-id	=> \@gb,
					-rettype	=> 'gb');

@hit_gis = split(m{\n},$factory->get_Response->content);

open (GB_TEXT,">>","hit_gb.txt") || die ("Error cannot write to hit_gb.txt: $!\n");

foreach $line (@hit_gis) {
	print GB_TEXT $line."\n";
}
close GB_TEXT;

open (IN2,"<","hit_gb.txt") || die ("Error cannot read from hit_gb.txt: $!\n");
@read = <IN2>;
close IN2;

open (OUT2,">>","hit_lineage.txt") || die ("Error cannot write to hit_lineage.txt:$!\n");

while ($read[$j]) {
	$line = $read[$j];

	chomp $line;
	if ($flag==0) {
		if ($line =~ /VERSION/) {
			$line =~ /VERSION\s+(\w+\.\d+)\s+GI:(\d+)/;
			$hit_gb = $1;
			$hit_gi = $2;
		}
		elsif ($line =~ /\s{2}ORGANISM/) {
			$line =~ s/  ORGANISM  //;
			@line = split(/ /,$line);
			$species = $line[1];
#			push(@species, $species);
			$flag=1;
		}
	}
	elsif ($flag==1) {
		if ($line !~ /REFERENCE/) {
				push(@organism, $line);
		}
		elsif ($line =~ /REFERENCE/) {
			$lineage = join(" ",@organism);
			$lineage =~ s/            //g;
			$lineage =~ s/.$/;/;
			$flag=0;
			print OUT2 "$hit_gi\t$hit_gb\t$lineage $species.\n";#test
			$hit_gb=();
			$hit_gi=();
			@organism=();
			$lineage=();
		}
	}
	$j++;
}
close OUT2;

open (IN3,"<","hit_lineage.txt") || die ("Error cannot read from hit_lienage.txt: $!\n");
@read2 = <IN3>;
close IN3;

open (OUT3,">>","hit_lineage.txt.mapped") || die ("Error cannot write to hit_lineage.txt.mapped: $!\n");

while ($read2[$k]) {
	$line = $read2[$k];
	chomp $line;

	@line = split(/\t/,$line);
	$hit_gb = $line[1];
#	print "hit_gb: $hit_gb\n";#test

	while ($gb[$l]) {
		$gb = $gb[$l];
		
		if ($gb =~ /$hit_gb/) {
			$query_map = $query_map[$l];
			print OUT3 "$query_map\t$line\n";	
		}
		$l++;
	}
	$l=0;
	$k++;
}


#!/usr/bin/perl
#March 31, 2011 by Terri Porter
#Script to use a list of full gb numbers to get gi numbers
#usage $ perl get_gi_for_each_gb.plx hit_gb.list

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

#declare array
my @in;
my @gb;
my @hit_gis;
my @organism;
my @read;

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@in= <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	push (@gb,$line);
	$i++;
}
#print "@gb\n";#test

$factory = Bio::DB::EUtilities-> new(	-eutil	=> 'efetch',
					-db	=>'nucleotide',
					-id	=> \@gb,
					-rettype	=> 'gb');

@hit_gis = split(m{\n},$factory->get_Response->content);

#my $scalar = scalar@hit_gis;
#print "scalar $scalar\n";#test
#$one_record = $hit_gis[1000];

open (GB_TEXT,">>","hit_gb.txt") || die ("Error cannot write to hit_gb.txt: $!\n");

foreach $line (@hit_gis) {
	print GB_TEXT $line."\n";
}
close GB_TEXT;
#print "@hit_gis\n";

#parse text formatted file and grab lineage, gi and gb

open (IN2,"<","hit_gb.txt") || die ("Error cannot read from hit_gb.txt: $!\n");
@read = <IN2>;
close IN2;

open (OUT2,">>","hit_lineage.txt") || die ("Error cannot write to hit_lineage.txt:$!\n");

while ($read[$j]) {
	$line = $read[$j];
#	print "$line\n";#test
	chomp $line;
	if ($flag==0) {
		if ($line =~ /VERSION/) {
			$line =~ /VERSION\s+(\w+\.\d+)\s+GI:(\d+)/;
			$hit_gb = $1;
			$hit_gi = $2;
			#print "got version\n";#test
			#print "hit_gb $hit_gb hit_gi $hit_gi\n";#test
		}
		elsif ($line =~ /\s{2}ORGANISM/) {
#			push(@organism, $line);
			$flag=1;
#			print "got organism\n";#test
		}
	}
	elsif ($flag==1) {
		if ($line !~ /REFERENCE/) {
				push(@organism, $line);
		}
		elsif ($line =~ /REFERENCE/) {
			$lineage = join(" ",@organism);
			$lineage =~ s/            //g;
			$flag=0;
			print OUT2 "$hit_gi\t$hit_gb\t$lineage\n";#test
			$hit_gb=();
			$hit_gi=();
			@organism=();
			$lineage=();
		}
	}
	$j++;
}
close OUT2;

#!/usr/bin/perl
#March 31, 2011 by Terri Porter
#Script to take gi\tgb\tlineage\n file and create gi_original\trank\n file
#usage $perl create_rank_summary.plx hit_lineage.txt gb_gi.query ITS.blast.parsed.all

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $original_gb;
my $hit_full_gb;
my $j=0;
my $k=0;
my $lineage;
my $rank;
my $l=0;
my $genus_maybe;
my $genus;
my $family;
my $order;
my $class;
my $phylum;
my $kingdom;
my $hit_gi;
my $hit_gb;
my $index;
my $filename;
my $next=0;
my $gi_part;
my $z=0;
my $end=0;
my $m=0;
my $map_gi;

#declare array
my @hit_lineage;
my @blast;
my @original_gb;
my @hit_full_gb;
my @line;
my @lineage;
my @lin;
my @gi_part;
my @map;
my @mapped_gis;

open (HIT_LINEAGE,"<",$ARGV[0]) || die ("Error cannot read from hit_lineage.txt: $!\n");
@hit_lineage = <HIT_LINEAGE>;
close HIT_LINEAGE;

open (MAP,"<",$ARGV[1]) || die ("Error cannot read from gb_gi.query: $!\n");
@map = <MAP>;
close MAP;

open (BLAST,"<",$ARGV[2]) || die ("Error cannot read from ITS.blast.parsed.all: $!\n");
@blast = <BLAST>;
close BLAST;

#keep track of original gb and corresponding hit_gb

while ($blast[$i]) {
	$line = $blast[$i];
	chomp $line;

	if ($line =~ /^\w+\s+\w+\|\w+\.\d+\|/) {
		$line =~ /^(\w+)\s+\w+\|(\w+\.\d+)\|/;
		$original_gb = $1;
		$hit_full_gb = $2;
#		print "original_gb $original_gb hit_full_gb $hit_full_gb\n";#test
		push(@original_gb,$original_gb);
		push(@hit_full_gb,$hit_full_gb);

		while ($map[$m]) {
			$line = $map[$m];
			chomp $line;
			$line =~ s/^\s+//;
			if ($line =~ /$original_gb/) {
				@line = split(/\t/,$line);
				$map_gi = $line[1];
				$map_gi =~ s/GI://;
				push(@mapped_gis,$map_gi);
			}
			$m++;
		}
		@line=();
		$m=0;
	}
	$i++;
}
#my $temp1 = scalar(@original_gb);
#my $temp2 = scalar(@hit_full_gb);
#print "original gb $temp1 hit full gb $temp2\n";#test
#print "map array @mapped_gis\n";#test
#parse lineage into appropriate ranks only

open (TEMP,">>","hit_lineage.parsed") || die ("Error cannot write to hit_lineage.temp: $!\n");

while ($hit_lineage[$l]) {
	$line = $hit_lineage[$l];
	chomp $line;
#	print "$line\n";#test

	@line = split(/;/,$line);
	
	my $temp7 = $line[0];
#	print "line index0 $l $temp7\n";#test

	while ($line[$z]) {
		$rank = $line[$z];
		$rank =~ s/^\s+//;
#		$rank =~ s/\.//;
		$rank =~ s/incertae sedis//;
#		$rank =~ s/mitosporic //;
#		$rank =~ s/unclassified //;

#		print $rank."\n";#test
		if ($next==1) {
			$rank =~ s/\.//;
			$genus_maybe = $rank;
			if ($genus_maybe =~ /^[A-Z]/) {
				$genus = $genus_maybe;
			}
			$next=0;
		}
		elsif ($rank =~ /aceae$/) {
			$family = $rank;
			$next=1;
		}
		elsif ($rank =~ /ales$/) {
			$order = $rank;
		}
		elsif ($rank =~ /Ichthyophonida/) {
			$order = $rank;
		}
		elsif ($rank =~ /mycetes$/) {
			$class = $rank;
		}
		elsif ($rank =~ /Ichthyosporea/) {
			$class = $rank;
		}
		elsif ($rank =~ /mycota$/) {
			$phylum = $rank;
		}
		elsif ($rank =~ /Fungi$/) {
			$kingdom = $rank;
		}
		elsif ($rank =~ /\t/) {
			$gi_part = $rank;
		#	$gi_part =~ /\d+\s+(\w+\.\d+)\s+(\w+)/;
			@gi_part = split(/\s+/,$gi_part);
		#	$hit_gi = $gi_part[0];
			$hit_gb = $gi_part[1];
		#	$kingdom = $2;
		#	print "$l hit_gb $hit_gb kingdom $kingdom\n";
		#	print "$l gi_part $gi_part\n";#test
		}
		$z++;
	}
	$z=0;
	print TEMP "$hit_gb\t$genus\t$family\t$order\t$class\t$phylum\t$kingdom\n";
	$l++;
	$gi_part=();
	$hit_gb=();
	$genus=();
	$family=();
	$order=();
	$class=();
	$phylum=();
	$kingdom=();
	$rank=0;
	@gi_part=();
}
close TEMP;

my $temp3 = $hit_lineage[17];
#my $temp4 = $line[17];
#print "hit lineage index17 $temp3\n";#test

print "Please enter what rank to summarize (ex. genus, family, order, class, phylum, kingdom):\n";
$rank = <STDIN>;
chomp $rank;

if ($rank eq 'genus') {
	$index = 1;
}
elsif ($rank eq 'family') {
	$index=2;
}
elsif ($rank eq 'order') {
	$index=3;
}
elsif ($rank eq 'class') {
	$index=4;
}
elsif ($rank eq 'phylum') {
	$index=5;
}
elsif ($rank eq 'kingdom') {
	$index=6;
}

open (LIN,"<","hit_lineage.parsed") || die ("Error cannot read from hit_lineage.temp: $!\n");
@lin = <LIN>;
close LIN;

$filename = "hit.".$rank;
open (RANK,">>",$filename) || die ("Error cannot write to $filename: $!\n");

while ($hit_full_gb[$j]) {
	$hit_full_gb = $hit_full_gb[$j];
#	print RANK "$j\t";#test

	while ($lin[$k]) {
		$line = $lin[$k];
		chomp $line;
		if ($end==0) {
			if ($line =~ /$hit_full_gb/) {
				@line = split(/\t/,$line);
				$rank = $line[$index];
				$original_gb = $mapped_gis[$j];
#				print RANK "$k\t";#test
				print RANK "$original_gb\t$rank\n";
				$end = 1;
			}
		}
		$k++;
	}
	if ($end==0) {
		print "$j\t$hit_full_gb\n";#test
	}
	$end=0;
	$k=0;
	$rank=();
	$original_gb=();
	$j++;
}
close RANK;
my $temp5 = $hit_full_gb[17];
my $temp6 = $lin[17];
#print "hit full gb index17 $temp5 lin index17 $temp6\n";#test

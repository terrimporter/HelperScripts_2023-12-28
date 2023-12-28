#!/usr/bin/perl
#May 28, 2011 by Terri Porter
#Work with clone_summary_best_rank.txt files from parse_clone_html_files_rank3.plx scripts
#usage recreate_hit_lineage.txt.mapped2.plx query_lineage.txt best_species.txt best_genus.txt best_family.txt best_order.txt best_class.txt best_phylum.txt best_kingdom.txt

use strict;
use warnings;

#declare var
my $i=0;
my $query_line;
my $map;
my $a=0;
my $kingdom_line;
my $kingdom;
my $flag=0;
my $b=0;
my $phylum_line;
my $phylum;
my $c=0;
my $class_line;
my $class;
my $d=0;
my $order_line;
my $order;
my $e=0;
my $family_line;
my $family;
my $f=0;
my $genus_line;
my $genus;
my $g=0;
my $species_line;
my $species;
my $species_stat;
my $genus_stat;
my $family_stat;
my $order_stat;
my $class_stat;
my $phylum_stat;
my $kingdom_stat;
my $stat;

#declare array
my @query;
my @species;
my @genus;
my @family;
my @order;
my @class;
my @phylum;
my @kingdom;
my @query_line;
my @kingdom_line;
my @phylum_line;
my @class_line;
my @order_line;
my @family_line;
my @genus_line;
my @species_line;

open (QUERY,"<",$ARGV[0]) || die ("Error cannot read query_lineage.txt: $!\n");
@query = <QUERY>;
close QUERY;

open (SPECIES,"<",$ARGV[1]) || die ("Error cannot read best_species.txt: $!\n");
@species = <SPECIES>;
close SPECIES;

open (GENUS,"<",$ARGV[2]) || die ("Error cannot read best_genus.txt: $!\n");
@genus = <GENUS>;
close GENUS;

open (FAMILY,"<",$ARGV[3]) || die ("Error cannot read best_family: $!\n");
@family = <FAMILY>;
close FAMILY;

open (ORDER,"<",$ARGV[4]) || die ("Error cannot read best_order: $!\n");
@order = <ORDER>;
close ORDER;

open (CLASS,"<",$ARGV[5]) || die ("Error cannot read best_class: $!\n");
@class = <CLASS>;
close CLASS;

open (PHYLUM,"<",$ARGV[6]) || die ("Error cannot read best_phylum: $!\n");
@phylum = <PHYLUM>;
close PHYLUM;

open (KINGDOM,"<",$ARGV[7]) || die ("Error cannot read best_kingdom: $!\n");
@kingdom = <KINGDOM>;
close KINGDOM;

open (OUT,">>","sap_hit_lineage.txt.mapped") || die ("Error cannot write to outfile: $!\n");

print "Please enter statistical cutoff (ex. 0-100):\n";
$stat = <STDIN>;
chomp $stat;

while($query[$i]) {
	$query_line = $query[$i];
	chomp $query_line;
	@query_line = split(/\t/,$query_line);
	$map = $query_line[0];
	print OUT "$map\t111\t111.1\tEukaryota;";
	while ($kingdom[$a]) {
		$kingdom_line = $kingdom[$a];
		chomp $kingdom_line;
		if ($kingdom_line =~ /$map/) {
			@kingdom_line = split(/\t/,$kingdom_line);
			$kingdom = $kingdom_line[1];
			$kingdom_stat = $kingdom_line[2];
			if ($kingdom_stat>=$stat) {
				print OUT " $kingdom; ";
				$flag=1;
			}
			else {
				print OUT "nil; ";
				$flag=1;
			}
		}
		$a++;
	}
	$a=0;
	@kingdom_line=();
	check_flag();

	while($phylum[$b]) {
		$phylum_line = $phylum[$b];
		chomp $phylum_line;
		
		if ($phylum_line =~ $map) {
			@phylum_line = split(/\t/,$phylum_line);
			$phylum = $phylum_line[1];
			$phylum_stat = $phylum_line[2];
			if ($phylum_stat>=$stat) {
				print OUT "$phylum; ";
				$flag=1;
			}
			else {
				print OUT "nil; ";
				$flag=1;
			}
		}
		$b++;
	}
	$b=0;
	@phylum_line=();
	check_flag();

	while ($class[$c]) {
		$class_line = $class[$c];
		chomp $class_line;
	
		if ($class_line =~ /$map/) {
			@class_line = split(/\t/,$class_line);
			$class = $class_line[1];
			$class_stat = $class_line[2];
			if ($class_stat>=$stat) {
				print OUT "$class; ";
				$flag=1;
			}
			else {
				print OUT "nil; ";
				$flag=1;
			}
		}
		$c++;
	}
	$c=0;
	@class_line=();
	check_flag();

	while ($order[$d]) {
		$order_line = $order[$d];
		chomp $order_line;

		if ($order_line =~ /$map/) {
			@order_line = split(/\t/,$order_line);
			$order = $order_line[1];
			$order_stat = $order_line[2];
			if ($order_stat>=$stat) {
				print OUT "$order; ";
				$flag=1;
			}
			else {
				print OUT "nil; ";
				$flag=1;
			}
		}
		$d++;
	}
	$d=0;
	@order_line=();
	check_flag();

	while ($family[$e]) {
		$family_line = $family[$e];
		chomp $family_line;

		if ($family_line =~ /$map/) {
			@family_line = split(/\t/,$family_line);
			$family = $family_line[1];
			$family_stat = $family_line[2];
			if ($family_stat>=$stat) {
				print OUT "$family; ";
				$flag=1;
			}
			else {
				print OUT "nil; ";
				$flag=1;
			}
		}
		$e++;
	}
	$e=0;
	@family_line=();
	check_flag();

	while($genus[$f]) {
		$genus_line = $genus[$f];
		chomp $genus_line;

		if ($genus_line =~ /$map/) {
			@genus_line = split(/\t/,$genus_line);
			$genus = $genus_line[1];
			$genus_stat = $genus_line[2];
			if ($genus_stat>=$stat) {
				print OUT "$genus;. ";
				$flag=1;
			}
			else {
				print OUT "nil;. ";
				$flag=1;
			}
		}
		$f++;
	}
	$f=0;
	@genus_line=();
	check_flag2();

	while($species[$g]) {
		$species_line = $species[$g];
		chomp $species_line;

		if ($species_line =~ /$map/) {
			@species_line = split(/\t/,$species_line);
			$species = $species_line[2];
			$species_stat = $species_line[3];
			if ($species_stat>=$stat) {
				print OUT "$species..";
				$flag=1;
			}
			else {
				print OUT "nil..";
				$flag=1;
			}
		}
		$g++;
	}
	$g=0;
	@species_line=();
	check_flag();
	print OUT "\n";

	$i++;
}

####################
sub check_flag {
	if ($flag==0) {
		print OUT "nil; ";
	}
	else {
		$flag=0;
	}
}
####################
sub check_flag2 {
	if ($flag==0) {
		print OUT"nil;. ";
	}
	else {
		$flag=0;
	}
}

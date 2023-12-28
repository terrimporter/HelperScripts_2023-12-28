#!/usr/bin/perl
#April 7, 2011 by Terri Porter
#Script to compare the hits from the 5' and 3' ends to see if the correct hits are the same or different
#usage $perl compare_hits.plx file.blastn.parsed.all hit_organism.filtered name.query

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $id_line;
my $hit_genus;
my $hit_species;
my $gi;
my $gb;
my $j=0;
my $line2;
my $ref_genus;
my $ref_species;
my $k=0;
my $line3;
my $hit_gb;
my $flag=0;
my $nothing = ();

#declare array
my @hit;
my @ref;
my @line;
my @id_line;
my @line2;
my @parsed;
my @line3;

open (HIT,"<",$ARGV[0]) || die ("Error cannot read file.blastn.parsed.all: $!\n");
@hit = <HIT>;
close HIT;

open (PARSED,"<",$ARGV[1]) || die ("Error cannot read hit_organism.filtered: $!\n");
@parsed = <PARSED>;
close PARSED;

open (REF,"<",$ARGV[2]) || die ("Error cannot read name.query: $!\n");
@ref = <REF>;
close REF;

open (OUT1,">>","genus.compared") || die ("Error cannot write to genus.compared: $!\n");
open (OUT2,">>","species.compared") || die ("Error cannot wrote to species.compared: $!|n");	

while ($hit[$i]) {
	$line = $hit[$i];
	chomp $line;
	
	if ($line =~ /^\S+\t\w+\|\w+\.\d+/) {
		$line =~ /^(\S+)\t\w+\|(\w+)\.\d+/;
		$gb = $1;
		$hit_gb = $2;

		while ($parsed[$k]) {
			$line3 = $parsed[$k];
			chomp $line3;
			if ($flag==0) {
			if ($line3 =~ /$hit_gb/) {
				@line3 = split(/\t/,$line3);
				$hit_genus = $line3[1];
				$hit_species = $line3[2];

				while($ref[$j]) {
					$line2 = $ref[$j];
					chomp $line2;
					
					if ($line2 =~ /$gb/) {
						@line2 = split(/\t/,$line2);
						$ref_genus = $line2[1];
						$ref_species = $line2[2];

						if ($ref_genus eq $hit_genus) {
							print OUT1 "$gb\t$ref_genus\t1\n";
						}
						elsif ($ref_genus ne $hit_genus) {
							print OUT1 "$gb\t$ref_genus\t0\n";
						}
						if ($ref_species eq $hit_species) {
							print OUT2 "$gb\t$ref_genus\t$ref_species\t1\n";
						}
						elsif ($ref_species ne $hit_species) {
							print OUT2 "$gb\t$ref_genus\t$ref_species\t0\n";
						}
						$flag=1;
					}
					
					$j++;
					@line2=();
					$ref_genus=();
					$ref_species=();
				}
				$j=0;
				$line2=();
			}
			}
			@line3=();
			$hit_genus=();
			$hit_species=();
			$k++;
		}
		$line3=();
		$k=0;
	}
	$flag=0;
	$i++;
	$line=();
	$gb=();
	$hit_gb=();
}
close OUT1;
close OUT2;


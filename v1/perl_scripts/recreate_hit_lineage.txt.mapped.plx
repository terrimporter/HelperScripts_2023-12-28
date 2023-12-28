#!/usr/bin/perl
#May 24, 2010 by Terri Porter
#Script to use MEGAN readid_taxonid.txt and output from R script using Bioconductor package genomes and fuction taxid2names to recreate the hit_lineage.txt.mapped format that can be parsed by create_rank_summary_4.plxs
#usage $perl recreate_hit_lineage.txt.mapped.plx megean_readid_taxid.txt output.txt

use strict;
use warnings;

#declare var
my $map;
my $map_id;
my $tax_id;
my $j=0;
my $i=0;
my $lineage;
my $organism;
my $line;
my $num_words;
my $flag=0;
my $genus_maybe;
my $species_maybe;
my $species;

#declare array
my @map;
my @rout;
my @map_line;
my @organism;

open (MAP,"<",$ARGV[0]) || die ("Error cannot read megan_readid_taxid.txt:$!\n");
@map = <MAP>;
close MAP;

open (ROUT,"<",$ARGV[1]) || die ("Error cannot read otuput.txt: $!\n");
@rout = <ROUT>;
close ROUT;

open (MEGAN,">>","megan_hit_lineage.txt.mapped") || die ("Error cannot write megan_hit_lineage.txt.mapped: $!\n");

while ($map[$i]) {
	$map = $map[$i];
	chomp $map;
#	print $map."\n";

	@map_line = split(/,/,$map);
	$map_id = $map_line[0];
	$tax_id = $map_line[1];
	$tax_id =~ s/^\s{1}//;

	while ($rout[$j]) {
		$line = $rout[$j];
		chomp $line;
#		print "$line\n";
		if ($flag==0) {
			if ($line =~ /$tax_id/) {
				$flag=1;
#				print "found line\n";
				$line =~ /"\d+"\s{1}$tax_id\s{1}"(.+)"\s{1}"(.+)"/g;
				$organism = $1;
				$lineage = $2;
				@organism = split(/ /,$organism);
				$num_words = scalar(@organism);
				$genus_maybe = $organism[0];
				$species_maybe = $organism[1];

				if ($num_words == 1) {
					print MEGAN "$map_id\t111\t111.1\t$lineage.\n";
				}

				elsif ($num_words ==2) {
					if ($genus_maybe =~ /^[A-Z]/) {
						if ($species_maybe =~ /^[a-z]/) {
							$species = $species_maybe;
							print MEGAN "$map_id\t111\t111.1\t$lineage;. $species..\n";
						}
					}
				}
			}		
		}
		$j++;
#		print $j."\n";
	}
	$flag=0;
	$line=();
	$organism=();
	$lineage=();
	$j=0;
	$i++;
	@map_line=();
	$map_id=();
	$tax_id=();
}
close MEGAN;

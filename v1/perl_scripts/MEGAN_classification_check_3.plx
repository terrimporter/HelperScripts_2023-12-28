#!/usr/bin/perl
#May 20, 2011 by Terri Porter
#Script to use query_lineage.txt with megan_readid_leaf.txt to check classifications
#usage MEGAN_classification_check_3.plx megan_readid_leaf.txt query_lineage.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $map;
my $gi;
my $gb;
my $leaf;
my $notassigned_counter=0;
my $words;
my $FungiMetazoa_counter=0;
my $species;
my $genus;
my $species_counter=0;
my $genus_counter=0;
my $family_counter=0;
my $order_counter=0;
my $class_counter=0;
my $phylum_counter=0;
my $kingdom_counter=0;
my $species_wrong=0;
my $genus_wrong=0;
my $family_wrong=0;
my $order_wrong=0;
my $class_wrong=0;
my $phylum_wrong=0;
my $kingdom_wrong=0;

#declare array
my @megan;
my @query;
my @line;
my @map;
my @leaf;


open (MEGAN,"<",$ARGV[0]) || die ("Error cannot read megan infile: $!\n");
@megan = <MEGAN>;
close MEGAN;

open (QUERY,"<",$ARGV[1]) || die ("Error cannot read query infile: $!\n");
@query = <QUERY>;
close QUERY;

while ($megan[$i]) {
	$line = $megan[$i];
	chomp $line;
	@line = split(/,/,$line);
	$map = $line[0];
	@map = split(/\|/,$map);
	$gi = $map[0];
	$gb = $map[1];
	$leaf = $line[1];
	$leaf =~ s/^\s{1}//;
	if ($leaf eq 'No hits') {
		$notassigned_counter++;
	}
	elsif ($leaf eq 'Not assigned') {
		$notassigned_counter++;
	}
#	elsif ($leaf =~ /incertae sedis/) {
#		$leaf =~ s/incertae sedis//;
#	}

	@leaf = split(/ /,$leaf);
	$words = scalar(@leaf);
	#genus level or higher
	if ($words == 1) {
		check_query();
	}
	#species level OR Fungi/Metazoa group
	elsif ($words == 2) {
		if ($leaf =~ /Fungi\/Metazoa group/) {
			$FungiMetazoa_counter++;
		}
		else {
			$species = $leaf[1];
			$genus = $leaf[0];
			check_query();
		}
	}
	elsif ($words == 3) { #genus sp. strain# OR class incertae sedis
		if ($leaf =~ /incertae sedis/) {
			$leaf = $leaf[0];
			check_query();
		}
		else {
			$genus = $leaf[0];
			check_query();
		}
	}
	elsif ($words == 4) {#genus cf. species strain# OR genus sp. strainletter strain#
		$genus = $leaf[0];
		check_query();
	}
	$i++;
	@line=();
	@map=();
	$gi=();
	$gb=();
	$leaf=();
	@leaf=();
	$words=();
	$genus=();
	$species=();
}

print "rank = correct\twrong\nspecies = $species_counter\t$species_wrong\ngenus = $genus_counter\t$genus_wrong\nfamily = $family_counter\t$family_wrong\norder = $order_counter\t$order_wrong\nclass = $class_counter\t$class_wrong\nphylum = $phylum_counter\t$phylum_wrong\nkingdom = $kingdom_counter\t$kingdom_wrong\nnot assigned = $notassigned_counter\n";

####################

sub check_query {

#declare var
my $j=0;
my $line2;
my $lineage;
my $num_lineage;
my $index;
my $rank_name;
my $higher_rank_name;
my $higher_rank_index;

#declare array
my @line2;
my @lineage;

	while ($query[$j]) {
		$line2 = $query[$j];
		chomp $line2;

		if ($line2 =~ /^$gi/) {
			if ($species) {
				if ($line2 =~ /$species\./) {
					$species_counter++;
					$genus=();
				}
				else {
					$species_wrong++;
				}
			}
			elsif ($genus) {
				if ($line2 =~ /$genus;\s{1}\w+\./) {
					$genus_counter++;
				}
				else {
					$genus_wrong++;
				}
			}
			else {
				if ($line2 =~ /$leaf/) {#track correct assignments
					if ($leaf =~ /Fungi/) {
						$kingdom_counter++;
					}
					elsif ($leaf =~ /mycota$/) {
						$phylum_counter++;
					}
					elsif ($leaf =~ /(mycetes$|Ichthyosporea)/) {
						$class_counter++;
					}
					elsif ($leaf =~ /(ales$|Ichthyophonida)/) {
						$order_counter++;
					}
					elsif ($leaf =~ /aceae$/) {
						$family_counter++;
					}
					elsif ($line2 =~ /$genus;\s{1}\w+\./) {
						$genus_counter++;
					}
					else {
						@line2 = split(/\t/,$line2);
						$lineage = $line2[2];
						@lineage = split(/;/,$lineage);
						$num_lineage = scalar(@lineage);
						$index = $num_lineage-1;

						while ($lineage[$index]) {
							if ($index>=0) {
								$rank_name = $lineage[$index];

								if ($rank_name =~ /$leaf/) {
									$higher_rank_index = $index-1;
									$higher_rank_name = $lineage[$higher_rank_index];
									if ($higher_rank_name =~ /(Fungi$|mycota$|ales$|aceae$)/) {
										$leaf = $higher_rank_name;
										$j--;#now recheck for a rank I'm interested in
									}
								}
								else {
									print "not counted $leaf\n";
								}
							}
							$index--;
						}

					}
				}
				else {#keep track of incorrect assignments to calculate coverage
					if ($leaf =~ /Viridiplantae\|Metazoa/) {
						$kingdom_wrong++;
					}
					elsif ($leaf =~ /(mycota$|Chlorophtyta|Streptophyta)/) {
						$phylum_wrong++;
					}
					elsif ($leaf =~ /(mycetes$|Ichthyosporea|opsida$|ophyceae$)/) {
						$class_wrong++;
					}
					elsif ($leaf =~ /(ales$|Ichthyophonida)/) {
						$order_wrong++;
					}
					elsif ($leaf =~ /aceae$/) {
						$family_wrong++;
					}
					else {
						$genus_wrong++;
					}
				}
			}
		}
		$j++;
	}
	$j=0;

}




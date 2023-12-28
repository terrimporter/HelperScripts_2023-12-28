#!/usr/bin/perl
# July 14/20 reformat output to make it more qiime format to accomodate the same taxon names at different ranks
#Apr. 11, 2018 retain species rank for CO1 v3git training set; handle 'superkingdom' rank and create 'cellularOrganisms'
#Jan. 5, 2017 edit to retain species rank for CO1v3 training set
#Oct. 24, 2016 add superkingdom to hash
#Aug. 11, 2016 add Protura
#March 22, 2013 by Terri Porter
#Script to create taxonomy file for Ribosomal Database Project Naive Bayesian Classifier
#usage perl make_NBC_taxonomy.plx taxid.parsed.awk.uniq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $cellularOrganisms;
my $superkingdom;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $termcounter=0;
my $j;
my $key;
my $value;
my $previous;
my $termcounter_previous;

#declare array
my @in;
my @value;
my @line;
my @previous;

#declare hash
my %lineage;
my %sorted;

open (IN, "<", $ARGV[0]) || die "Error cannot open taxid.parsed.awk.uniq: $!\n";
@in = <IN>;
close IN;

#hash lineage table and create relationships
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);
	
	$cellularOrganisms = $line[0];
	$cellularOrganisms =~ s/\s+//g; #removing trailing space
	
	$superkingdom = $line[1];
	$superkingdom =~ s/\s+//g; #remove any spaces
	
	$kingdom = $line[2];
	$kingdom =~ s/\s+//g;
	if ($kingdom eq 'undef') {
		$kingdom = $kingdom.'_'.$superkingdom;
	}
	$phylum = $line[3];
	$phylum =~ s/\s+//g;
	if ($phylum eq 'undef') {
		$phylum = $phylum.'_'.$kingdom;
	}
#	if ($phylum eq 'Acanthocephala') { # phylum of worms, genus of insects
#		$phylum = $phylum.'_'.$kingdom;
#	}
#	if ($phylum eq 'Ctenophora') { ### phylum of Metazoa and genus of arthropods
#		$phylum = $phylum.'_'.$kingdom;
#	}
#	if ($phylum eq 'Actinobacteria') { ### phylum and class of Bacteria
#		$phylum = $phylum.'_'.$kingdom;
#	}
#	if ($phylum eq 'Gemmatimonadetes') { ### phylum and class of Bacteria
#		$phylum = $phylum.'_'.$kingdom;
#	}
#	if ($phylum eq 'Deferribacteres') { ### phylum and class of Bacteria
#		$phylum = $phylum.'_'.$kingdom;
#	}

	$class = $line[4];
	$class =~ s/\s+//g;
	if ($class eq 'undef') {
		$class = $class.'_'.$phylum;
	}
#	if ($class eq 'Deferribacteres') { # phylum and class of Bacteria
#		$class = $class.'_'.$phylum;
#	}


	$order = $line[5];
	$order =~ s/\s+//g;
	if ($order eq 'undef') {
		$order = $order.'_'.$class;
	}
#	if ($order eq 'Plecoptera') { ### Plecoptera is genus of moths and order of stoneflies
#		$order = $order.'_'.$class;
#	}
#	if ($order eq 'Protura') { ### Protura is class and order
#		$order = $order.'_'.$class;
#	}
#	if ($order eq 'Diplura') { ### Diplura is class and order
#		$order = $order.'_'.$class;
#	}
	if ($order eq 'Pristiformes/Rhiniformesgroup') { ### fix / and missing space
		$order = 'Pristiformes_Rhiniformes_group';
	}
#	if ($order eq 'Parachela') { ### order of tardigrades and genus of fish
#		$order = $order.'_'.$class;
#	}
	if ($order eq 'BacteroidetesOrderII.Incertaesedis') { ### fix spacing
		$order = 'Bacteroidetes_Order_II_Incertae_sedis';
	}
#	if ($order eq 'Pygophora') { # genus flies, order crustaceans
#		$order = $order.'_'.$class;
#	}

	$family = $line[6];
	$family =~ s/\s+//g;
	if ($family eq 'undef') {
		$family = $family.'_'.$order;
	}
	if ($family eq 'Cepheidae') { ### family of mites and jellyfish
		$family = $family.'_'.$order;
	}
	if ($family eq 'Chilodontidae') { ### family of fish and molluscs
		$family = $family.'_'.$order;
	}
	$genus = $line[7];
	$genus =~ s/\s+//g;

	if ($genus eq 'undef') {
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Caviria') { # genus in Lymantriidae in NCBI, Erebidae in BOLD
		$family = 'Lymantriidae';
	}
	if ($genus eq 'Antachara') { # genus in Noctuidae in NCBI, Erebidae in BOLD
		$family = 'Noctuidae';
	}
	if ($genus eq 'Eisenia') { # genus brown algae and segmented worms
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Lactarius') { # genus bony fishes & basidiomycetes
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Ricinus') { # genus lice and eudicots
		$genus = $genus.'_'.$family;
	}
#	if ($genus eq 'Pygophora') { # genus flies, order crustaceans
#		$genus = $genus.'_'.$family;
#	}
	if ($genus eq 'Eutrapela') { ### Eutrapela is found in family Geometridae and Tenebrionidae
		$genus = $genus.'_'.$family;
	}
#	if ($genus eq 'Plecoptera') {### Plecoptera is genus of moths and order of stoneflies
#		$genus = $genus.'_'.$family;
#	}
	if ($genus eq 'Alaria') {### Alaria is genus of stramenopiles and platyhelminths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Automolus') {### Automolus is a genus of birds and beetles
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Achlya') {### Achlya is a genus of oomycetes and lepidoptera
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Ozophora') {### Ozophora is a genus of insects and red algae
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Acrotylus') {### Acrotylus is a genus of red algae and insects
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Lessonia') {### Lessonia is a genus of birds and brown algae
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Chondracanthus') {### Chondracanthus is a genus of red algae and copepods
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Elachista') { ### genus of stramenopiles and moths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Bostrychia') { ### genus of stoneflies and birds
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Roya') { ### genus of green algae and molluscs
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Rhizophagus') { ### genus of beetle and fungus
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Nemastoma') { ### genus of arthropod and red alga
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Ariopsis') { ### genus of plants and catfish
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Pontogeneia') { ### genus of fungi and amphipods
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Chondria') { ### genus of red algae and beetles
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Ptilophora') { ### genus of red algae and moths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Karenia') { ### genus of alveolata and cicada
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Lobophora') { ### genus of stramenopiles and moths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Candida') { ### 2 clades of Candida yeasts
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Dacrydium') { ### genus of conifer and bivalve
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Trigonostomum') { ### genus of flatworms and beetles
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Nectria') { ### genus of fungi and worms
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Dracunculus') { ### genus of plants and nematodes
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Ganonema') { ### genus of red algae and caddisflies
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Bouchetia') { ### genus of plant and sea snail
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Agardhiella') { ### genus of red algae and land snails
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Solieria') { ### genus of red algae and flies
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Diplura') { ### genus of arthropods and stramenopiles
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Zonaria') { ### genus of stramenopiles and gastropods
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Vestia') { ### genus of plants and snails
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Hillia') { ### genus of plants and moths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Olea') { ### genus of plants and sea slugs
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Grania') { ### genus of annelid worms and red algae
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Atractomorpha') { ### genus of plants and insects
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Bartramia') { ### genus of mosses and birds
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Heterococcus') { ### genus of stramenopiles and arthropods
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Tephrosia') { ### genus of plants and moths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Rondeletia') { ### genus of plants and fish
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Turbinaria') { ### genus of cnidarians and stramenopiles
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Mastophora') { ### genus of arthropods and red algae
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Drymonia') { ### genus of plants and moths
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Planococcus') { ### genus of bacteria and insects
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Paracoccus') { ### genus of bacteria and mealy bugs
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Rhodococcus') { ### genus of bacteria and scale insects
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Bacillus') { ### genus of bacteria and stick insects
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Darwinella') { ### genus of beetles and sponges
		$genus = $genus.'_'.$family;
	}

	$species = $line[8];
	$species =~ s/^\s+//g; #remove any preceeding whitespace
	$species =~ s/\s+/_/g; #change spaces to underscores just in case
	$species =~ s/_$//g; #if underscore at end of line, remove it
	if ($species eq 'undef') {
		$i++;
		next;
	}
	elsif ($species =~ /Pelophila_borealis/) {
		$phylum = 'Arthropoda';
		$class = 'Insecta';
		$order = 'Coleoptera';
		$family = 'Carabidae';
		$genus = 'Pelophila';

	}
	elsif ($species =~ /coxendix$/) {
		$species = 'Apallates_coxendix'; # weird character/spacing
	}

	if (exists $lineage{$cellularOrganisms}) {
		
		if (exists $lineage{$superkingdom}) {
	
			if (exists $lineage{$kingdom}) {

				if (exists $lineage{$phylum}) {
		
					if (exists $lineage{$class}) {
			
						if (exists $lineage{$order}) {
					
							if (exists $lineage{$family}) {
						
								if (exists $lineage{$genus}) {

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
							else {
								create_family();

								if (exists $lineage{$genus}) {
									$i++;
									next;
								}
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
						else {
							create_order();

							if (exists $lineage{$family}) {
								$i++;
								next;
							}
							else {
								create_family();

								if (exists $lineage{$genus}) {
									$i++;
									next;
								}
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
					}
					else {
						create_class();

						if (exists $lineage{$order}) {
							$i++;
							next;
						}
						else {
							create_order();
	
							if (exists $lineage{$family}) {
								$i++;
								next;
							}
							else {
								create_family();
	
								if (exists $lineage{$genus}) {
									$i++;
									next;
								}
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
					}
				}
				else {
					create_phylum();

					if (exists $lineage{$class}) {
						$i++;
						next;
					}
					else {
						create_class();

						if (exists $lineage{$order}) {
							$i++;
							next;
						}
						else {
							create_order();

							if (exists $lineage{$family}) {
								$i++;
								next;
							}
							else {
								create_family();

								if (exists $lineage{$genus}) {
									$i++;
									next;
								}
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
					}
				}
			}
			else {
				create_kingdom();

				if (exists $lineage{$phylum}) {
					$i++;
					next;
				}				
				else {
					create_phylum();

					if (exists $lineage{$class}) {
						$i++;						
						next;
					}
					else {
						create_class();

						if (exists $lineage{$order}) {
							$i++;
							next;
						}
						else {
							create_order();

							if (exists $lineage{$family}) {
								$i++;
								next;							
							}
							else {
								create_family();

								if (exists $lineage{$genus}) {
									$i++;
									next;
								}								
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
					}
				}
			}
		}
		else {
			create_superkingdom();

			if (exists $lineage{$kingdom}) {
				$i++;
				next;
			}
			else {
				create_kingdom();

				if (exists $lineage{$phylum}) {
					$i++;
					next;
				}				
				else {
					create_phylum();

					if (exists $lineage{$class}) {
						$i++;						
						next;
					}
					else {
						create_class();

						if (exists $lineage{$order}) {
							$i++;
							next;
						}
						else {
							create_order();

							if (exists $lineage{$family}) {
								$i++;
								next;							
							}
							else {
								create_family();

								if (exists $lineage{$genus}) {
									$i++;
									next;
								}								
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {
		$termcounter++;
		$j=$termcounter-1;
		$lineage{$cellularOrganisms} = $termcounter."*r__".$cellularOrganisms."*".$j."*0*cellularOrganisms";
		
		if (exists $lineage{$superkingdom}) {
			$i++;
			next;
		}
		else {
			create_superkingdom();

			if (exists $lineage{$kingdom}) {
				$i++;
				next;
			}
			else {
				create_kingdom();

				if (exists $lineage{$phylum}) {
					$i++;
					next;
				}
				else {
					create_phylum();

					if (exists $lineage{$class}) {
						$i++;
						next;
					}
					else {
						create_class();

						if (exists $lineage{$order}) {
							$i++;
							next;
						}
						else {
							create_order();

							if (exists $lineage{$family}) {
								$i++;
								next;
							}
							else {
								create_family();

								if (exists $lineage{$genus}) {
									$i++;
									next;
								}
								else {
									create_genus();

									if (exists $lineage{$species}) {
										$i++;
										next;
									}
									else {
										create_species();
									}
								}
							}
						}
					}
				}
			}
		}
	}
	$i++;
	$line=();
	@line=();
	$cellularOrganisms=();
	$superkingdom=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
}
$i=0;
$termcounter=();

#sort the hash by termcounter into a new hash indexed by termcounter
while ( ($key, $value) = each (%lineage) ) {
	@value = split(/\*/, $value);
	$termcounter = $value[0];
	$sorted{$termcounter} = $value;
}
$key=();
$value=();

#sort keys then print

open (OUT, ">>", "testNBC.taxonomy") || die "Error cannot open outfile:$!\n";

foreach $key (sort {$a <=> $b} keys %sorted) {
	$value = $sorted{$key};
	print OUT $value."\n";
}
close OUT;

##########

sub create_species {

	$termcounter++;
	$previous = $lineage{$genus};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$species} = $termcounter."*s__".$species."*".$termcounter_previous."*8*species";

}

##########

sub create_genus {

	$termcounter++;
	$previous = $lineage{$family};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$genus} = $termcounter."*g__".$genus."*".$termcounter_previous."*7*genus";

}

###########

sub create_family {

	$termcounter++;
	$previous = $lineage{$order};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$family} = $termcounter."*f__".$family."*".$termcounter_previous."*6*family";

}

############

sub create_order {

	$termcounter++;
	$previous = $lineage{$class};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$order} = $termcounter."*o__".$order."*".$termcounter_previous."*5*order";

}

#############

sub create_class {

	$termcounter++;
	$previous = $lineage{$phylum};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$class} = $termcounter."*c__".$class."*".$termcounter_previous."*4*class";

}

##############

sub create_phylum {

	$termcounter++;
	$previous = $lineage{$kingdom};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$phylum} = $termcounter."*p__".$phylum."*".$termcounter_previous."*3*phylum";

}

###############

sub create_kingdom {

	$termcounter++;
	$previous = $lineage{$superkingdom};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$kingdom} = $termcounter."*k__".$kingdom."*".$termcounter_previous."*2*kingdom";

}

###############

sub create_superkingdom {

	$termcounter++;
	$previous = $lineage{$cellularOrganisms};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$superkingdom} = $termcounter."*sk__".$superkingdom."*".$termcounter_previous."*1*superkingdom";

}

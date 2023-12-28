#!/usr/bin/perl

# May 9, 2019 by Teresita M. Porter
# Script to parse BOLD_DR.fasta (formatted for RDP classifier), only keep records with a complete species field
# get rid of last BIN field for now
# only keep fully identified species because they tend to have filled out lineages, ignore BIN's for now
# create a list of id_poorqual.txt
# create list of gb_seq
# create list of gb_taxid
# create filtered_taxids.txt
# create taxid.parsed ## should have undef when BOLD taxon not in NCBI
# usage perl filter_BOLD_species.plx BOLD_DR.fasta

use strict;
use warnings;

# declare var
my $id_poorqual_out = "id_poorqual.txt";
my $gb_seq_out = "gb_seq.map";
my $gb_taxid_out = "gb_taxid.map";
my $filtered_taxids = "filtered_taxids.txt";
my $taxid_parsed = "taxid.parsed";

my $namespath = "/home/terri/ncbi-blast-2.9.0+/db/names.dmp";
my $nodespath = "/home/terri/ncbi-blast-2.9.0+/db/nodes.dmp";

my $i=0;
my $line;
my $species;
my $words;
my $j;
my $seq;
my $speciesEpithet; #from BOLD
my $id;
my $id_root;
my $length;
my $root;
my $superkingdom;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $count;
my $taxidlist;
my $original_taxid;
my $parent;
my $rank;
my $taxid;
my $name;
my $nameclass;
my $oldtaxid;
my $newtaxid;
my $genus_taxid;
my $family_taxid;
my $order_taxid;
my $class_taxid;
my $phylum_taxid;
my $kingdom_taxid;
my $superkingdom_taxid;
my $genusChoice;
my $taxidChoice;

# declare array
my @in;
my @header;
my @species;
my @words;
my @id_root;
my @names;
my @nodes;
my @merged;
my @taxids;
my @line;

# declare hash
my %taxid_name; 	#key=taxid, 	value=name, 		should be unique
my %name_taxid; 	#key=name, 		value=taxid, 		names could have more than one taxid
my %taxid_parent; 	#key=taxid, 	value=parent_taxid
my %taxid_rank; 	#key=taxid, 	value=rank
my %old_new;		#key=oldtaxid, 	value=newtaxid		old taxids are sometimes merged into newer ones, check here if taxid not found

open (NAMES, "<", $namespath) || die "Cannot open names.dmp: $!\n";
@names = <NAMES>;
close NAMES;

# hash names for quick lookups
while ($names[$i]) {
	$line = $names[$i];
	chomp $line;

	@line = split(/\t\|\t/,$line);
	$taxid = $line[0];
	$name = $line[1];
	$nameclass = $line[3];
	$nameclass =~ s/\t\|//g; #remove final pipe

	if ($nameclass eq "scientific name") { # in the future also consider synonyms ##########
		$taxid_name{$taxid} = $name; # should always be unique

		if (exists $name_taxid{$name}){ # name could have more than one taxid, do not overwrite
			$original_taxid = $name_taxid{$name};
			$taxidlist = $original_taxid.";".$taxid;
			$name_taxid{$name} = $taxidlist;
		}
		else {
			$name_taxid{$name} = $taxid;
		}
	}
	$i++;
}
$i=0;

#print $name_taxid{"Bruchomorpha tristis"}."\n"; #test

open (NODES, "<", $nodespath) || die "Cannot open nodes.dmp: $!\n";
@nodes = <NODES>;
close NODES;

# hash nodes for quick lookups
while ($nodes[$i]) {
	$line = $nodes[$i];
	chomp $line;

	@line = split(/\t\|\t/,$line);
	$taxid = $line[0];
	$parent = $line[1];
	$rank = $line[2];

	$taxid_parent{$taxid} = $parent;
	$taxid_rank{$taxid} = $rank;

	$i++;
}
$i=0;

# read in the cleaned up BOLD seqs in RDP FASTA format
open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

# create a bunch of new outfiles for make_NBC_fasta_cellularOrganisms_species_omitpoorqual.plx
open (ID, ">>", $id_poorqual_out) || die "Error cannot open id_poorqual_out: $!\n";
open (GB_SEQ, ">>", $gb_seq_out) || die "Error cannot open gb_seq_out: $!\n";
open (GB_TAXID, ">>", $gb_taxid_out) || die "Error cannot open gb_taxid_out: $!\n";
open (FILT, ">>", $filtered_taxids) || die "Error cannot open filtered_taxids:$!\n";
open (PARS, ">>", $taxid_parsed) || die "Error cannot open taxid_parsed:$!\n";

# start parsing the infile
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@header = split(";",$line);
		$id_root = shift(@header); # shift first element off top of array
		@id_root = split(/ /,$id_root); # BOLD sampldIDs can have spaces!?!
		$root = pop(@id_root); #pop last element off end of array 'cellularOrganisms'
		$id = join ' ',@id_root; #join whats left in the sampleID (otherwise GenBank accession or BOLD BIN would be fine)
		$id =~ s/^>//g; # remove leading greater than sign from FASTA
		$id =~ s/ /_/g; #replace spaces with underscore so taxid.parsed not messed up
		$superkingdom = shift(@header);
		$kingdom = shift(@header);
		$phylum = shift(@header);
		$class = shift(@header);
		$order = shift(@header);
		$family = shift(@header);
		$genus = shift(@header);
		$species = shift(@header);

		##### Handle weird exceptions #####
		if (length $species) {
			if ($species =~ m/^Erigorgus stenotus/) { # NCBI uses Erigorgus as subgenus but Gravenhorstia as genus, BOLD uses Erigorgus as genus
				$genus = "Gravenhorstia";
			}
			if ($species =~ m/^Murchieonia/) {
				$species =~ s/Murchieonia/Murchieona/;
#				print "species:".$species."\n"; #test
			}
			if ($species =~ m/^Cyclura lechriodes/) {
				$species =~ s/Cyclura lechriodes/Amphitorna lechriodes/; # Cyclura is a synonym for Amphitorna, but it's also a genus of iguana
				$genus = "Amphitorna";
			}
			if ($species =~ m/^Urodon suturalis/) {
				$species =~ s/Urodon suturalis/Bruchela suturalis/; # Urodon is a synonym for Bruchela, but it's also a genus of plants
				$genus = "Bruchela";
			}
			if ($species =~ m/^Urodon conformis/) {
				$species =~ s/Urodon conformis/Bruchela conformis/; # Urodon is a synonym for Bruchela, but it's also a genus of plants
				$genus = "Bruchela";
			}
			if ($species =~ m/^Urodon rufipes/) {
				$species =~ s/Urodon rufipes/Bruchela rufipes/; # Urodon is a synonym for Bruchela, but it's also a genus of plants
				$genus = "Bruchela";
			}

		}
		if (length $genus) {
			if ($genus eq 'Murchieonia') { # BOLD data release has Murchieonia but it's Murchieona everywhere else.
				$genus = "Murchieona";
			}
		}
		if (length $class) {
			if ($class eq 'Actinopterygii') { # BOLD calls Actinopterygii as a class but NCBI thinks its a superclass, use Actinopteri to be consistent with NCBI
				$class = "Actinopteri";
			}
		}
		if (length $class) {
			if ($class eq 'Elasmobranchii') { # BOLD calls Elasmobranchii a class but NCBI thinks its a subclass, NCBI uses Chondrichthyes
				$class = "Chondrichthyes";
			}
		}

 		#if there is a BOLD BIN in this field, means there isn't a proper species epithet
		if (length $species && $species =~ /BOLD/) {
			print ID "$id\n"; #start a list of ids not to process
			$i+=2;
			next;
		}
		elsif (length $species) { #make sure species field isn't empty
#			print "species found: $species\n";
			@species = split(/ /,$species);
			$words = scalar(@species);
			$speciesEpithet = $species[1];

			if ($words > 2) {
				print ID "$id\n"; #don't process records with species that contain more than 2 words
					$i+=2;
					next;
			}
			elsif ($words == 2 && $speciesEpithet =~ /^[A-Z]/) {
				print ID "$id\n"; #species epithet should never begin with a capital else probably not real species ID
				$i+=2;
				next;
			}
			elsif ($species =~ /(sp\d+|sp\.|aff\.|cf\.|nr\.)/) {
				print ID "$id\n"; #don't process insufficiently identified records (for now)
				$i+=2;
				next;
			}
			elsif ($line =~ /\;\;/) {
				print ID "$id\n"; # skip any records with missing or shifted taxonomy (problem parsing tsv)
				$i+=2;
				next;
			}
			else {
				$j = $i+1;
				$seq = $in[$j];
				chomp $seq;

				print GB_SEQ "$id\t$seq\n";

				if (exists $name_taxid{$species}) {
					$taxidlist = $name_taxid{$species};
					check_for_dups($taxidlist, "species"); #check for multiple taxids in list, pick just one at the right rank, print to outfiles
#					print "species still? $species\n"; #test
					# reconcile BOLD taxonomy against NCBI taxonomy
					#keep genus and species as-is whether or not it's in Genbank

					if (exists $name_taxid{$genus}) {
						$genus_taxid = $name_taxid{$genus}; # possible for there to be more than one taxid for a name

						if (length $genus_taxid) {
							$genus_taxid = check_for_dups($genus_taxid, "genus");

							if (length $genus_taxid) {
								$count=0; #set counter to prevent deep recursions
								$family = rank_check("family",$genus_taxid); # get ncbi taxon
							}
							else {
								print "Couldn't find parent taxid for genus $genus\n";
							}
						}
						else {
							$family = "undef";
						}
						$family_taxid = $name_taxid{$family};

						if (length $family_taxid) {
							$family_taxid = check_for_dups($family_taxid, "family");
							$count=0;

							if (length $family_taxid) {
								$order = rank_check("order",$family_taxid); # get ncbi taxon
							}
							else {
								$family="undef";
							}
						}
						else {
							$family="undef";
						}
						$order_taxid = $name_taxid{$order};
						if (length $order_taxid) {
							$order_taxid = check_for_dups($order_taxid, "order");
							$count=0;
							if (length $order_taxid) {
								$class = rank_check("class",$order_taxid); # get ncbi taxon
							}
							else {
								$order="undef";
							}
						}
						else {
							$order="undef";
						}
						$class_taxid = $name_taxid{$class};
						if (length $class_taxid) {
							$class_taxid = check_for_dups($class_taxid, "class");
							$count=0;
							if (length $class_taxid) {
								$phylum = rank_check("phylum",$class_taxid); # get ncbi taxon
							}
							else {
								$class="undef";
							}
						}	
						else {
							$class="undef";
						}
						$phylum_taxid = $name_taxid{$phylum};
						if (length $phylum_taxid) {
							$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
							$count=0;
							if (length $phylum_taxid) {
								$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
							}
							else {
								$phylum="undef";
							}
						}
						else {
							$phylum="undef";
						}
						$kingdom_taxid = $name_taxid{$kingdom};
						if (length $kingdom_taxid) {
							$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
							$count=0;
							if (length $kingdom_taxid) {
								$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
							}
							else {
								$kingdom="undef";
							}
						}
						else {
							$kingdom="undef";
						}
						print PARS $taxidlist."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
						print FILT $taxidlist."\n";
						$i+=2;
						next;
					}
					else {
						#genus doesn't exist in NCBI!
#						$genus="";
						$family="undef"; 
						# proceed with order...
						if (exists $name_taxid{$order}) { #assumd BOLD order is good and grab remaining lineage
							$order_taxid = $name_taxid{$order};
							$order_taxid = check_for_dups($order_taxid, "order");
							$count=0;
							if (length $order_taxid) {
								$class = rank_check("class",$order_taxid); # get ncbi taxon
							}
							else {
								$order="undef";
							}
							if (exists $name_taxid{$class}) {
								$class_taxid = $name_taxid{$class};
								$class_taxid = check_for_dups($class_taxid, "class");
								$count=0;
								if (length $class_taxid) {
									$phylum = rank_check("phylum",$class_taxid); # get ncbi taxon
								}
								else {
									$class="undef";
								}
							}
							else {
								$class="undef";
							}
							if (exists $name_taxid{$phylum}) {
								$phylum_taxid = $name_taxid{$phylum};
								$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
								$count=0;
								if (length $phylum_taxid) {
									$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
								}
								else {
									$phylum = "undef";
								}
							}
							else {
								$phylum="undef";
							}
							if (exists $name_taxid{$kingdom}) {
								$kingdom_taxid = $name_taxid{$kingdom};
								$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
								$count=0;
								if (length $kingdom_taxid) {
									$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
								}
								else {
									$kingdom="undef";
								}
							}
							else {
								$kingdom="undef";
							}
							print PARS $taxidlist."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
							print FILT $taxidlist."\n";
							$i+=2;
							next;
						}
						else {
							$order = "undef";
							if (exists $name_taxid{$class}) {
								$class_taxid = $name_taxid{$class};
								$class_taxid = check_for_dups($class_taxid, "class");
								$count=0;
								if (length $class_taxid) {
									$phylum = rank_check("phylum",$class_taxid); # get ncbi taxon
								}
								else {
									$class = "undef";
								}
								if (exists $name_taxid{$phylum}) {
									$phylum_taxid = $name_taxid{$phylum};
									$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
									$count=0;
									if (length $phylum_taxid) {
										$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
									}
									else {
										$phylum="undef";
									}
								}
								else {
									$phylum="undef";
								}
								if (exists $name_taxid{$kingdom}) {
									$kingdom_taxid = $name_taxid{$kingdom};
									$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
									$count=0;
									if (length $kingdom_taxid) {
										$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
									}
									else {
										$kingdom="undef";
									}
								}
								else {
									$kingdom="undef";
								}
								print PARS $taxidlist."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
								print FILT $taxidlist."\n";
								$i+=2;
								next;
							}
							else {
								$class = "undef";
								if (exists $name_taxid{$phylum}) {
									$phylum_taxid = $name_taxid{$phylum};
									$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
									$count=0;
									if (length $phylum_taxid) {
										$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
									}
									else {
										$phylum="undef";
									}
									if (exists $name_taxid{$kingdom}){
										$kingdom_taxid = $name_taxid{$kingdom};
										$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
										$count=0;
										if (length $kingdom_taxid) {
											$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
										}
										else {
											$kingdom="undef";
										}
									}
									else {
										$kingdom="undef";
									}
									print PARS $taxidlist."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
									print FILT $taxidlist."\n";
									$i+=2;
									next;
								}
								else {
									$phylum="undef";
									if (exists $name_taxid{$kingdom}) {
										$kingdom_taxid = $name_taxid{$kingdom};
										$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
										$count=0;
										if (length $kingdom_taxid) {
											$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
										}
										else {
											$kingdom="undef";
										}
										print PARS $taxidlist."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
										print FILT $taxidlist."\n";
										$i+=2;
										next;
									}
									else {
										$kingdom = "undef";
										unless (exists $name_taxid{$superkingdom}) {
											$superkingdom="undef";
											print PARS $taxidlist."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
											print FILT $taxidlist."\n";
											$i+=2;
										}
									}
								}
							}
						}
					}
				}
				else {
					# create a placeholder taxid for BOLD species not in GenBank
					$taxid = "TEMP".$i;
					# add to hash so it can be looked up in future entries
					$name_taxid{$species} = $taxid;
#					print "species3 $species\n"; #test
					print GB_TAXID "$id\t$taxid\n";
					print FILT "$taxid\n";

					if (exists $name_taxid{$genus}) {
						$genus_taxid = $name_taxid{$genus};
						$genus_taxid = check_for_dups($genus_taxid, "genus");
						$count=0; #set counter to prevent deep recursions
						if (length $genus_taxid) {
							$family = rank_check("family",$genus_taxid); # get ncbi taxon
						}
						else {
							print "Can't find taxid for genus $genus\n";
						}
						if (exists $name_taxid{$family} ) {
							$family_taxid = $name_taxid{$family};
							$family_taxid = check_for_dups($family_taxid, "family");
							$count=0;
							if (length $family_taxid) {
								$order = rank_check("order",$family_taxid); # get ncbi taxon
							}
							else {
								$family="undef";
							}
						}
						else {
							$family = "undef";
						}
						if (exists $name_taxid{$order}) {
							$order_taxid = $name_taxid{$order};
							$order_taxid = check_for_dups($order_taxid, "order");
							$count=0;
							if (length $order_taxid) {
								$class = rank_check("class",$order_taxid); # get ncbi taxon
							}
							else {
								$order="undef";
							}
						}
						else {
							$order = "undef";
						}
						if (exists $name_taxid{$class}) {
							$class_taxid = $name_taxid{$class};
							$class_taxid = check_for_dups($class_taxid, "class");
							$count=0;
							if (length $class_taxid) {
								$phylum = rank_check("phylum",$class_taxid); # get ncbi taxon
							}
							else {
								$class="undef";
							}
						}
						else {
							$class="undef";
						}
						if (exists $name_taxid{$phylum}) {
							$phylum_taxid = $name_taxid{$phylum};
							$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
							$count=0;
							if (length $phylum_taxid) {
								$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
							}
							else {
								$phylum="undef";
							}
						}
						else {
							$phylum="undef";
						}
						if (exists $name_taxid{$kingdom}) {
							$kingdom_taxid = $name_taxid{$kingdom};
							$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
							$count=0;
							if (length $kingdom_taxid) {
								$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
							}
							else {
								$kingdom="undef";
							}
						}
						else {
							$kingdom="undef";
						}
						print PARS $taxid."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
						print FILT $taxid."\n";
						$i+=2;
						next;
					}
					else {
						#handle this, family="", keep class+
						$family = "undef"; # BOLD often classifies genera in families different from NCBI so keep as undef
						if (exists $name_taxid{$order}) { #assumd BOLD order is good and grab remaining lineage
							$order_taxid = $name_taxid{$order};
							$order_taxid = check_for_dups($order_taxid, "order");
							$count=0;
							if (length $order_taxid) {
								$class = rank_check("class",$order_taxid); # get ncbi taxon
							}
							else {
								$order="undef";
							}
							if (exists $name_taxid{$class}) {
								$class_taxid = $name_taxid{$class};
								$class_taxid = check_for_dups($class_taxid, "class");
								$count=0;
								if (length $class_taxid) {
									$phylum = rank_check("phylum",$class_taxid); # get ncbi taxon
								}
								else {
									$class="undef";
								}
							}
							else {
								$class="undef";
							}
							if (exists $name_taxid{$phylum}) {
								$phylum_taxid = $name_taxid{$phylum};
								$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
								$count=0;
								if (length $phylum_taxid) {
									$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
								}
								else {
									$phylum="undef";
								}
							}
							else {
								$phylum="undef";
							}
							if (exists $name_taxid{$kingdom}) {
								$kingdom_taxid = $name_taxid{$kingdom};
								$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
								$count=0;
								if (length $kingdom_taxid) {
									$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
								}
								else {
									$kingdom="undef";
								}
							}
							else {
								$kingdom="undef";
							}
							print PARS $taxid."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
							print FILT $taxid."\n";
							$i+=2;
							next;
						}
						else {
							$order = "undef";
							if (exists $name_taxid{$class}) {
								$class_taxid = $name_taxid{$class};
								$class_taxid = check_for_dups($class_taxid, "class");
								$count=0;
								if (length $class_taxid) {
									$phylum = rank_check("phylum",$class_taxid); # get ncbi taxon
								}
								else {
									$class="undef";
								}
								if (exists $name_taxid{$phylum}) {
									$phylum_taxid = $name_taxid{$phylum};
									$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
									$count=0;
									if (length $phylum_taxid) {
										$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
									}
									else {
										$phylum="undef";
									}
								}
								else {
									$phylum="undef";
								}
								if (exists $name_taxid{$kingdom}) {
									$kingdom_taxid = $name_taxid{$kingdom};
									$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
									$count=0;
									if (length $kingdom_taxid) {
										$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
									}
									else {
										$kingdom="undef";
									}
								}
								else {
									$kingdom="undef";
								}
								print PARS $taxid."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
								print FILT $taxid."\n";
								$i+=2;
								next;
							}
							else {
								$class = "undef";
								if (exists $name_taxid{$phylum}) {
									$phylum_taxid = $name_taxid{$phylum};
									$phylum_taxid = check_for_dups($phylum_taxid, "phylum");
									$count=0;
									if (length $phylum_taxid) {
										$kingdom = rank_check("kingdom",$phylum_taxid); # get ncbi taxon
									}
									else {
										$phylum="undef";
									}
									if (exists $name_taxid{$kingdom}){
										$kingdom_taxid = $name_taxid{$kingdom};
										$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
										$count=0;
										if (length $kingdom_taxid) {
											$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
										}
										else {
											$kingdom="undef";
										}
									}
									else {
										$kingdom="undef";
									}
									print PARS $taxid."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
									print FILT $taxid."\n";
									$i+=2;
									next;
								}
								else {
									$phylum="undef";
									if (exists $name_taxid{$kingdom}) {
										$kingdom_taxid = $name_taxid{$kingdom};
										$kingdom_taxid = check_for_dups($kingdom_taxid, "kingdom");
										$count=0;
										if (length $kingdom_taxid) {
											$superkingdom = rank_check("superkingdom",$kingdom_taxid); # get ncbi taxon
										}
										else {
											$kingdom="undef";
										}
										print PARS $taxid."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
										print FILT $taxid."\n";
										$i+=2;
										next;
									}
									else {
										$kingdom = "undef";
										unless (exists $name_taxid{$superkingdom}) {
											$superkingdom="undef";
											print PARS $taxid."\t".$root."\t".$superkingdom."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
											print FILT $taxid."\n";
											$i+=2;
											next;
										}
									}
								}
							}
						}
					}
					$taxid=();
				}
			}
		}
		else {
			print ID "$id\n"; #don't process insufficiently identified records
			$i+=2;
			next;
		}
	}
	else {
		$i++;
		next;
	}
	$genusChoice=();
	$taxidChoice=();
	$genus_taxid=();
	$family_taxid=();
	$order_taxid=();
	$class_taxid=();
	$phylum_taxid=();
	$kingdom_taxid=();
	$taxid=();

}
$i=0;
$taxidChoice=();
close ID;
close GB_SEQ;
close GB_TAXID;
close FILT;
close PARS;

#############################

sub rank_check {

# subroutine to recursively crawl up taxonomy tree until desired rank is found and return taxon name
# no more than 35 levels permitted before "" returned

my $targetRank = shift;
my $taxidLocal = shift;

my $parent;
my $rank;
my $value;

	if ($count == 0) {

		if (length $taxidLocal) {

			if (exists $taxid_parent{$taxidLocal}) {
				$parent = $taxid_parent{$taxidLocal}; # check taxid up one level

				if (exists $taxid_rank{$parent}) {
					$rank = $taxid_rank{$parent};

					if ($rank eq $targetRank) {
						$value = $taxid_name{$parent};
						return $value;
					}
					else {
						$count++;
						$value = recursive_rank_check($targetRank, $parent);
						return $value;
					}
				}
				else {
					print "can't find parent $parent in taxid_rank\n";
					$value = "";
					return $value;
				}
			}
			else {
				print "can't find taxidLocal $taxidLocal in taxid_parent\n";
				$value = "";
				return $value;
			}
		}
		else {
			print "taxidLocal $taxidLocal is empty\n";
			$value = "";
			return $value;
		}
	}
	else {
		$value = "";
		return $value;
		print "count wasn't 0\n";
	}
}

#############################

sub recursive_rank_check {

my $targetRank = shift;
my $taxidLocal = shift;

my $parent;
my $rank;
my $value;

	if ($count < 35) { # fish have a lot of levels in NCBI!
		$parent = $taxid_parent{$taxidLocal};
		$rank = $taxid_rank{$parent};

		if ($rank eq $targetRank){
			$value = $taxid_name{$parent};
			return $value;
		}
		else {
			$count++;
			recursive_rank_check($targetRank, $parent);
		}
	}
	else {
		$value = "";
		return $value;
	}
}

############################################################

sub check_for_dups {

#Subroutine to use a name to check for taxid(s), can be a single taxid or a list, check rank and phylum assignment for each taxid, if it's good then return the taxid (original or correct one from the list)

my $taxidLocal = shift;
my $targetRank = shift;

#print "taxidLocal $taxidLocal, targetRank $targetRank\n";

my $rankChoice;
my $phylumCheck;
my $lengthLocal;

	if (length $taxidLocal) {
		if ($taxidLocal =~ /;/) { #some names can be associated with more than one taxid (especially genus+)
			@taxids = split(/;/,$taxidLocal);
			$lengthLocal = scalar(@taxids);

			foreach $taxid (@taxids) {
				$rankChoice = $taxid_rank{$taxid};

				if ($targetRank eq $rankChoice) {

					# ensure this taxid comes from the correct phylum!
					if ($targetRank eq "superkingdom" | $targetRank eq "kingdom" | $targetRank eq "phylum") {
						$taxidChoice = $taxid;
					}
					else {
						$count=0;
						$phylumCheck = rank_check("phylum",$taxid);
	
						if ($phylum eq $phylumCheck) { #compare phylum in BOLD taxonomy with phylum from taxid chosen from taxidlist
							$taxidChoice = $taxid;
						}
						else {
							#wrong phylum, don't do anything here	
							next;
						}
					}
				}
				else {
					#wrong genus, don't do anything here
					next;
				}
			}
		
			if ($targetRank eq "species") { # no need to return a taxid at the species rank, should never get here
				if (length $taxidChoice) {
					print GB_TAXID "$id\t$taxidChoice\n";
					print FILT "$taxidChoice\n";
				}
				else {
					# couldn't find the right taxid, don't bother processing this entry
					print "Couldn't find species level taxid $taxidChoice for id $id\n"; 
					$taxidChoice="";
					return $taxidChoice;
				}
			}
			else {
				if (length $taxidChoice) {
					return $taxidChoice;
				}
				else {
					$taxidChoice="";
					return $taxidChoice;
				}
			}
		}
		else {
#			return $taxidLocal; #it's not a list, just one taxid

			if ($targetRank eq "species") { # no need to return a taxid at the species rank, but need to print to outfiles
				if (length $taxidLocal) {
					print GB_TAXID "$id\t$taxidLocal\n";
					print FILT "$taxidLocal\n";
#					print "taxidLocal $taxidLocal\n"; #test, never get here???
					return $taxidLocal;
				}
				else {
					# couldn't find the right taxid, don't bother processing this entry
					print "Couldn't find species level taxid2 $taxidLocal\n"; #should never really get here?!?
				}
			}
			else {
				if (length $taxidLocal) {
					return $taxidLocal;
				}
				else {
					$taxidLocal="";
					return $taxidLocal;
				}
			}
		}
	}
	else {
		print "taxid $taxid is empty2\n";
		$taxidLocal="";
		return $taxidLocal;
	}
}

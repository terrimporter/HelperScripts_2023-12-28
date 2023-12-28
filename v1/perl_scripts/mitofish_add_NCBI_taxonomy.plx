# Teresita M. Porter, March 17, 2020
# Script to add NCBI taxonomy to the MitoFish species
# USAGE perl mitofish_add_NCBI_taxonomy.plx 12S.fasta names.dmp nodes.dmp

use strict;
use warnings;
use Data::Dumper;

# var
my $i=0;
my $line;
my $species;
my $grepline;
my $taxid;
my $acc1;
my $acc2;
my $acc;
my $next_taxid;
my $current_rank;
my $name;
my $lineage="";
my $counter=0;
my $j;
my $seq;
my $rank;
my $sciname;
my $old_species;
my $outfile = "testNBC.fasta";
my $outfile2 = "taxid.parsed";
my $species_taxid;
my $newlineage;
my $value;

# array
my @fasta;
my @line;
my @grep;
my @grepline;
my @grepline2;
my @names;
my @nodes;
my @node_search;
my @name_search;
my @matches;
my @lineage;
my @taxids;

# hash
my %species_names; # key = species binomial, value = NCBI taxid
my %target_ranks = ("species", 1,
					"genus", 1,
					"family", 1,
					"order", 1,
					"class", 1,
					"phylum", 1,
					"kingdom", 1,
					"superkingdom", 1
					);  # key = ranks, value = 1
my %acc_species; # key = acc, value = species binomial
my %acc_seq; # key = acc, value = seq
my %species_lineage; # key = binomial species, value = lineage
my %names; # key = taxid, value = name
my %nodes_nexttaxid; # key = taxid, value = next_taxid
my %nodes_rank; #key = taxid, value = rank
my %taxid_species; # key = taxid, value = species
my %lineage; # key = rank, value = name
my %species_mapping; # key = badly formatted species, value = correctly formatted species
my %synonym; # key = species, value = taxid
my %taxidparsed; #key = lineage, value = 1

# hash names.dmp
open (IN, "<", $ARGV[1]) || die "Error cannot open names.dmp:$!\n";
@names = <IN>;
close IN;

while ($names[$i]) {
	$line = $names[$i];
	chomp $line;

	@line = split(/\|/, $line);
	$taxid = $line[0];
	$taxid =~ s/\s+//g;
	$name = $line[1];
	$name =~ s/^\s+//g;
	$name =~ s/\s+$//g;
	$sciname = $line[3];
	$sciname =~ s/^\s+//g;
	$sciname =~ s/\s+$//g;

	# only add if it's a scientific name (ignore synonyms)
	if ($sciname eq "scientific name") {
		$names{$taxid} = $name;
		$species_names{$name} = $taxid;
	}
	elsif ($sciname eq "synonym") {
		$synonym{$name} = $taxid;
	}
	elsif ($sciname eq "includes") {
		$synonym{$name} = $taxid;
	}
	else {
		$i++;
		next;
	}

	$i++;

}
$i=0;

print "Finished parsing names.dmp\n";

open (IN, "<", $ARGV[0]) || die "Error cannot open 12S.fasta infile: $!\n";
@fasta = <IN>;
close IN;

# grab species names from FASTA headers
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~/^>/) {
		$line =~ s/_{2}/_/g;
		$line =~ s/\s+$//g;
		@line = split(/_/,$line);
		$acc1 = shift @line; #remove NC
		$acc1 =~ s/^>//g;
		$acc2 = shift @line; # remove rest of accession
		$species = join " ", @line;
		$species =~ s/ sp / sp\. /;
		$species =~ s/ cf / cf\. /;
		$species =~ s/mirror/'mirror'/;
		$species =~ s/singuonensis/'singuonensis'/;
		$species =~ s/southern/'southern'/;
		$species =~ s/var$/var\./;
		# get accession to handle the really messed up ones
		$acc = $acc1."_".$acc2;
		if ($acc eq "NC_031807") {
			$species = "(Cyprinus carpio 'mirror' x Cyprinus carpio 'singuonensis') x Carassius auratus red var.";
		}
		elsif ($acc eq "NC_028198") {
			$species = "Coregonus sp. 'cluncaformis'";
		}
		elsif ($acc eq "NC_036033") {
			$species = "Megalobrama amblycephala x (Megalobrama amblycephala x Parabramis pekinensis)";
		}
		elsif ($acc eq "NC_013615") {
			$species = "Leptocephalus sp. 'type II larva'";
		}
		elsif ($acc eq "NC_036034") {
			$species = "Megalobrama amblycephala x (Megalobrama amblycephala x Megalobrama terminalis)";
		}
		elsif ($acc eq "NC_015142") {
			$species = "Carassius auratus ssp. 'Pingxiang'";
		}
		elsif ($acc eq "NC_026543") {
			$species = "Cyprinus carpio 'Furong' x Carassius auratus red var.";
		}

		# found if it's a scientific name
		if (exists $species_names{$species}) {
			$taxid = $species_names{$species};
#			print "Found species $species in the hash\n";
		}
		# species may be a 'synonym' or 'includes' and not in names hash 
		elsif (exists $synonym{$species}) {
			$old_species = $species;
			@grep = `grep "$species" $ARGV[1]`;
			$grepline = $grep[0];
			@grepline = split(/\|/, $grepline);
			$taxid = $grepline[0];
			$taxid =~ s/\s+$//g;
			@grep = `grep "^$taxid " $ARGV[1]`;
			foreach $grepline (@grep) {
				@grepline = split(/\|/, $grepline);
				$sciname = $grepline[3];
				$sciname =~ s/^\s+//g;
				$sciname =~ s/\s+$//g;
				if ($sciname eq "scientific name") {
					$species = $grepline[1];
					$species =~ s/^\s+//g;
					$species =~ s/\s+$//g;
					$species = $species_mapping{$old_species};
				}
			}
#			print "Found species $species in the synonym hash\n";
		}
		# species may appear in a list among subspecies that start with same 'species' pattern
		else {
			@grep = `grep "$species" $ARGV[1]`;
			foreach $grepline (@grep) {
				@grepline = split(/\|/, $grepline);
				$sciname = $grepline[3];
				$sciname =~ s/^\s+//g;
				$sciname =~ s/\s+$//g;
				if ($sciname eq "scientific name") {
					$species = $grepline[1];
					$species =~ s/^\s+//g;
					$species =~ s/\s+$//g;
					$taxid = $grepline[0];
					$taxid=~ s/\s+$//g;
				}
			}
#			print "Found the species from partial grep match\n";
		}

		if (!defined($taxid)) {
			print "Still can't find taxid for $acc $species \n";
		}
	
		$taxid_species{$taxid} = $species;
		$acc_species{$acc} = $species;
		$j = $i+1;
		$seq = $fasta[$j];
		chomp $seq;
		$acc_seq{$acc} = $seq;
	}
	$i++;

	@grep=();
	$taxid=();
	$species=();

}
$i=0;

print "Finished parsing FASTA\n";

# hash nodes.dmp
open (IN, "<", $ARGV[2]) || die "Error cannot open nodes.dmp: $!\n";
@nodes = <IN>;
close IN;

while ($nodes[$i]) {
	$line = $nodes[$i];
	chomp $line;

	@line = split(/\|/, $line);
	$taxid = $line[0];
	$taxid =~ s/\s+$//g;
	$next_taxid = $line[1];
	$next_taxid =~ s/\s+//g;
	$rank = $line[2];
	$rank =~ s/\s+//g;
	$nodes_nexttaxid{$taxid} = $next_taxid;
	$nodes_rank{$taxid} = $rank;

	$i++;

}
$i=0;

print "Finished parsing nodes.dmp\n";

# for each taxid, grab taxonomy from nodes.dmp
while (($taxid, $species) = each(%taxid_species)) {
	$species = $taxid_species{$taxid};
	$species_taxid = $taxid;

	($next_taxid, $current_rank) = node_search($species_taxid);
	
	if ($current_rank eq "species") {
		$lineage{"species"} = $species;
	}
	elsif ($taxid eq "86966") { # Anguilla nebulosa nebulosa
		$old_species = $species;
		$species = "Anguilla nebulosa";
#		delete($taxid_species{$taxid});
#		$taxid_species{"86965"} = $species;
		$next_taxid = 7935;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "90313") { # cherry salmon
		$old_species = $species;
		$species = "Oncorhynchus masou";
#		delete($taxid_species{$taxid});
#		$taxid_species{"8020"} = $species;
		$next_taxid = 8016;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "749190") { # Cyprinus carpio 'color'
		$old_species = $species;
		$species = "Cyprinus carpio";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7962"} = $species;
		$next_taxid = 7961;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "149978") { # Sarcocheilichthys variegatus microoculus
		$old_species = $species;
		$species = "Sarcocheilichthys variegatus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"149977"} = $species;
		$next_taxid = 149979;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "291359") { # Cheilopogon pinnatibarbatus japonicus
		$old_species = $species;
		$species = "Cheilopogon pinnatibarbatus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"307590"} = $species;
		$next_taxid = 293908;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	} 
	elsif ($taxid eq "334711") { # Gymnocypris przewalskii ganzihonensis
		$old_species = $species;
		$species = "Gymnocypris przewalskii";
#		delete($taxid_species{$taxid});
#		$taxid_species{"75348"} = $species;
		$next_taxid = 75347;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "227976") { # Salmo trutta trutta
		$old_species = $species;
		$species = "Salmo trutta";
#		delete($taxid_species{$taxid});
#		$taxid_species{"8032"} = $species;
		$next_taxid = 8028;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "630221") { # Cyprinus carpio carpio
		$old_species = $species;
		$species = "Cyprinus carpio";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7962"} = $species;
		$next_taxid = 7961;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "72053") { # Oncorhynchus clarkii henshawi
		$old_species = $species;
		$species = "Oncorhynchus clarkii";
#		delete($taxid_species{$taxid});
#		$taxid_species{"30962"} = $species;
		$next_taxid = 8016;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "61125") { # Anguilla bicolor bicolor
		$old_species = $species;
		$species = "Anguilla bicolor";
#		delete($taxid_species{$taxid});
#		$taxid_species{"61124"} = $species;
		$next_taxid = 7935;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "286774") { # Phoxinus oxycephalus jouyi
		$old_species = $species;
		$species = "Phoxinus oxycephalus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"109666"} = $species;
		$next_taxid = 42662;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "357291") { # Rhynchocypris percnurus mantschuricus
		$old_species = $species;
		$species = "Rhynchocypris percnurus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"219678"} = $species;
		$next_taxid = 933989;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "1045198") { # Cyprinus carpio 'xingguonensis'
		$old_species = $species;
		$species = "Cyprinus carpio";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7962"} = $species;
		$next_taxid = 7961;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "8242") { # Thunnus thynnus thynnus
		$old_species = $species;
		$species = "Thunnus thynnus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"8237"} = $species;
		$next_taxid = 8234;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "1088067") { # Polypterus palmas buettikoferi
		$old_species = $species;
		$species = "Polypterus palmas";
#		delete($taxid_species{$taxid});
#		$taxid_species{"645109"} = $species;
		$next_taxid = 8290;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "204997") { # Sander vitreus vitreus
		$old_species = $species;
		$species = "Sander vitreus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"283036"} = $species;
		$next_taxid = 283033;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "763460") { # Carassius auratus ssp. 'Pingxiang'
		$old_species = $species;
		$species = "Carassius auratus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7957"} = $species;
		$next_taxid = 7956;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "487156") { # Trichiurus lepturus nanhaiensis
		$old_species = $species;
		$species = "Trichiurus lepturus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"13733"} = $species;
		$next_taxid = 13732;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "489037") { # Micropterus salmoides salmoides
		$old_species = $species;
		$species = "Micropterus salmoides";
#		delete($taxid_species{$taxid});
#		$taxid_species{"27706"} = $species;
		$next_taxid = 27705;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "178879") { # Polypterus senegalus senegalus
		$old_species = $species;
		$species = "Polypterus senegalus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"55291"} = $species;
		$next_taxid = 8290;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "86961") { # Anguilla australis australis
		$old_species = $species;
		$species = "Anguilla australis";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7940"} = $species;
		$next_taxid = 7935;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "2570877") { # Cyprinus carpio x Megalobrama amblycephala
		$old_species = $species;
		$species = "Cyprinidae intergeneric hybrids";
#		delete($taxid_species{$taxid});
#		$taxid_species{"564289"} = $species;
		$next_taxid = 7953;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "80800") { # Rhodeus ocellatus kurumeus
		$old_species = $species;
		$species = "Rhodeus ocellatus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"80800"} = $species;
		$next_taxid = 58326;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "8021") { # Oncorhynchus masou ishikawae
		$old_species = $species;
		$species = "Oncorhynchus masou";
#		delete($taxid_species{$taxid});
#		$taxid_species{"8020"} = $species;
		$next_taxid = 8016;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "103832") { # Anguilla bengalensis labiata
		$old_species = $species;
		$species = "Anguilla bengalensis";
#		delete($taxid_species{$taxid});
#		$taxid_species{"103813"} = $species;
		$next_taxid = 7935;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "1045269") { # Cyprinus carpio haematopterus
		$old_species = $species;
		$species = "Cyprinus carpio";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7962"} = $species;
		$next_taxid = 7961;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "319226") { # Rhynchocypris percnurus sachalinensis
		$old_species = $species;
		$species = "Rhynchocypris percnurus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"219678"} = $species;
		$next_taxid = 933989;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "244512") { # Petrochromis trewavasae ephippium
		$old_species = $species;
		$species = "Petrochromis trewavasae";
#		delete($taxid_species{$taxid});
#		$taxid_species{"244518"} = $species;
		$next_taxid = 28817;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "341118") { # Cobitis melanoleuca granoei
		$old_species = $species;
		$species = "Cobitis melanoleuca";
#		delete($taxid_species{$taxid});
#		$taxid_species{"425482"} = $species;
		$next_taxid = 47718;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "412663") { # Oncorhynchus masou 'Biwa'
		$old_species = $species;
		$species = "Oncorhynchus masou";
#		delete($taxid_species{$taxid});
#		$taxid_species{"8020"} = $species;
		$next_taxid = 8016;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "341105") { # Pseudogastromyzon fasciatus jiulongjiangensis
		$old_species = $species;
		$species = "Pseudogastromyzon fasciatus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"425505"} = $species;
		$next_taxid = 241466;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "481698") { # Tanakia himantegus himantegus
		$old_species = $species;
		$species = "Tanakia himantegus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"481777"} = $species;
		$next_taxid = 80803;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "86962") { # Anguilla australis schmidti
		$old_species = $species;
		$species = "Anguilla australis";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7940"} = $species;
		$next_taxid = 7935;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ( $taxid eq "57445") { # Maylandia zebra complex is a 'species group'
		$old_species = $species;
		$species = "Maylandia zebra complex";
#		delete($taxid_species{$taxid});
#		$taxid_species{"57445"} = $species;
		$next_taxid = 143623;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "1115630") { # Brachymystax lenok tsinlingensis
		$old_species = $species;
		$species = "Brachymystax lenok";
#		delete($taxid_species{$taxid});
#		$taxid_species{"62067"} = $species;
		$next_taxid = 62066;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "749191") { # Cyprinus carpio 'wuyuanensis'
		$old_species = $species;
		$species = "Cyprinus carpio";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7962"} = $species;
		$next_taxid = 7961;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "175781") { # Glyptothorax fokiensis fokiensis
		$old_species = $species;
		$species = "Glyptothorax fokiensis";
#		delete($taxid_species{$taxid});
#		$taxid_species{"175780"} = $species;
		$next_taxid = 175779;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "145527") { # Carassius auratus auratus
		$old_species = $species;
		$species = "Carassius auratus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"7957"} = $species;
		$next_taxid = 7956;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "442960") { # Gadus chalcogrammus chalcogrammus
		$old_species = $species;
		$species = "Gadus chalcogrammus";
#		delete($taxid_species{$taxid});
#		$taxid_species{"1042646"} = $species;
		$next_taxid = 8048;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "61126") { # Anguilla bicolor pacifica
		$old_species = $species;
		$species = "Anguilla bicolor";
#		delete($taxid_species{$taxid});
#		$taxid_species{"61124"} = $species;
		$next_taxid = 7935;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	elsif ($taxid eq "173242") { # Oncorhynchus masou formosanus
		$old_species = $species;
		$species = "Oncorhynchus masou";
#		delete($taxid_species{$taxid});
#		$taxid_species{"8020"} = $species;
		$next_taxid = 8016;
		$lineage{"species"} = $species;
		$species_mapping{$old_species} = $species;
	}
	else {
		print "Warning, first taxid $taxid doesn't represent a species $current_rank: $!\n";
	}

	($name, $taxid) = name_search($next_taxid);
	($next_taxid, $current_rank) = node_search($taxid);

	if ($current_rank eq "genus") {
		$lineage{"genus"} = $name;
	}
	else {
		if (exists $target_ranks{$current_rank}) {
			$lineage{$current_rank} = $name;
		}
	}
	
	($name, $taxid) = name_search($next_taxid);
	($next_taxid, $current_rank) = node_search($taxid);


	# loop through the rest of the lineage up to 40 nodes deep
	while ($counter < 40) {
		if ($current_rank eq "superkingdom") {
			$lineage{"superkingdom"} = $name;
			$counter=40;
		}
		elsif (exists $target_ranks{$current_rank}) {
			$lineage{$current_rank} = $name;
			($name, $taxid) = name_search($next_taxid);
			($next_taxid, $current_rank) = node_search($taxid);
		}
		else {
			($name, $taxid) = name_search($next_taxid);
			($next_taxid, $current_rank) = node_search($taxid);
		}
		$counter++;
	}
	$counter=0;

	# account for undefined lineages
	unless (exists $lineage{"superkingdom"}) {
		$lineage{"superkingdom"} = "undef";
	}
	unless (exists $lineage{"kingdom"}) {
		$lineage{"kingdom"} = "undef_".$lineage{"superkingdom"};
	}
	unless (exists $lineage{"phylum"}) {
		$lineage{"phylum"} = "undef_".$lineage{"kingdom"};
	}
	unless (exists $lineage{"class"}) {
		$lineage{"class"} = "undef_".$lineage{"phylum"};
	}
	unless (exists $lineage{"order"}) {
		$lineage{"order"} = "undef_".$lineage{"class"};
	}
	unless (exists $lineage{"family"}) {
		$lineage{"family"} = "undef_".$lineage{"order"};
	}
	unless (exists $lineage{"genus"}) {
		$lineage{"genus"} = "undef_".$lineage{"family"};
	}
	unless (exists $lineage{"species"}) {
		print "Warning, no species given for species_taxid $species_taxid!\n";
	}

	$lineage = $lineage{"superkingdom"}.";".$lineage{"kingdom"}.";".$lineage{"phylum"}.";".$lineage{"class"}.";".$lineage{"order"}.";".$lineage{"family"}.";".$lineage{"genus"}.";".$lineage{"species"};

	# replace space with underscore
	$lineage =~ s/\s{1}/_/g;

	$species_lineage{$species} = $lineage;

	%lineage=();

}

#print Dumper(\%acc_species);

#=begin
open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";
open (OUT2, ">>", $outfile2) || die "Error cannot open outfile2:$!\n";

# for each acc, get species, for each species get lineage, for each acc get seq
# print testNBC.fasta
while (($acc, $old_species) = each(%acc_species)) {
	$old_species = $acc_species{$acc};

	# check if species name was badly formatted 
	if (exists $species_mapping{$old_species}) {
		$species = $species_mapping{$old_species};
	}
	else {
		$species = $old_species;
	}

	# get lineage for species
	if (exists $species_lineage{$species}) {
		$lineage = $species_lineage{$species};
	}
	else {
		print "Warning, couldn't find lineage for acc $acc old species $old_species\n";
	}

	# get seq for acc
	if (exists $acc_seq{$acc}) {
		$seq = $acc_seq{$acc};
	}
	else {
		print "Warning, couldn't find seq for acc $acc\n";
	}

	print OUT ">$acc\tcellularOrganisms;$lineage\n$seq\n";
	
	@lineage = split(/;/,$lineage);
	$newlineage = join "\t", @lineage;
	$newlineage = "cellularOrganisms\t$newlineage";

#	print OUT2 "cellularOrganisms\t$newlineage";
	$taxidparsed{$newlineage} = 1;

}
close OUT;

# print out unique lineages only
# print taxid.parsed
while (($newlineage,$value) = each(%taxidparsed)) {
	print OUT2 $newlineage."\n";
}
close OUT2;

##############################################

sub node_search {

$taxid = $_[0];

	# grab the ndoe line based on taxid
	if (exists $nodes_rank{$taxid}) {
		$current_rank = $nodes_rank{$taxid};
		if (exists $nodes_nexttaxid{$taxid}) {
			$next_taxid = $nodes_nexttaxid{$taxid};
		}
		else {
			print "Warning, cannot find $taxid in nodes_nexttaxid hash\n";
		}
	}
	else {
		print "Warning, cannot find $taxid in nodes_rank hash\n";
	}

	return ($next_taxid, $current_rank);

}

##############################################

sub name_search {

	$next_taxid = $_[0];

	# grab name line based on taxid
	if (exists $names{$next_taxid}) {
			$name = $names{$next_taxid};
	}
	else {
		print "Warning, cannot find $next_taxid in names hash\n";
	}
		
	$taxid = $next_taxid;
	return($name, $taxid);

}


#!/usr/bin/perl
# Jan. 18, 2020 account for 18S sequences from Chloroplast or Mitochondria amongst problem taxa
# Jan. 17/20 edit to include Bacterial outgroups
# Jan. 15/20 edit to accomodate a mapping file of problematic lineages in SILVA fasta file where a species can have more than one unique lineage, and missing genera in SILVA taxonomy file
#July 27, 2018 include problem_species.txt file (from check_for_SIVLA_inconsistencies.plx)
#July 25, 2018 by Teresita M. Porter
#Script to fix FASTA file so only whole taxon ranks are represented (domain, kingdom, phylum, class, order, family, genus, species)
#convert U's to T's, in lower case
#ensure output is in strict FASTA format (no line breaks in seq)
#USAGE perl fix_SILVA_fasta_headers.plx tax_slv_ssu_138.txt SILVA_138_SSURef_NR99_tax_silva.fasta.gz problem_map.txt silva.centroids.gz

use strict;
use warnings;
use Data::Dumper;

#var
my $i=0;
my $line;
my $path;
my $lineage;
my $abundance;
my $rank;
my $acc;
my $fulltaxon;
my $taxon;
my $prevtaxon="";
my $newlineage;
my $outfile = "testNBC.fasta";
my $seq="";
my $euk=0; #flag to grab just Eukaryota
my $domain;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $prob=0; #flag record to skip, 0 = ok, 1 = in %problem, 2 = using mapping from in %fixproblem
my $good_seq_counter=0;
my $salvaged_seq_counter=0;
my $bad_seq_counter=0;
my $bact_seq_counter=0;
my $originalLineage;

#arrays
my @tax;
my @fasta;
my @line;
my @lineage;
my @taxon;
my @problem;
my @centroids;

#hashes
my %map; #key=path, value=rank
my %duplicates; #key=taxon, value=1
my %observed; #key = observed rank in SILVA_132_SSUParc_tax_silva.fasta.gz, value = observed taxon name
my %problem; #outerkey = species, innerkey = lineage, value = abundance
my %fixproblem; #key = species, value = lienage

# handle duplicate taxon names (same name across different ranks and/or lineages)
# may represent duplicates because of SSU from mitochondria or chloroplasts in addition to the nucleus
# add to this list when errors pop up during RDP classifier training

%duplicates= (
		"Lepismatidae;" => 1,
		"Thelebolus;"  => 1, # tax issues
		"Syllides;" => 1, # GU357624 is an annelid not an arthropod
		"Cassiopea;" => 1, # tax issues
		"Lamproderma;" => 1, # tax issues
		"Lacrymaria;" => 1,
		"Rhizophydium;" => 1, # tax issues
		"Athelia;" => 1, # tax issues
		"Rana;" => 1, # EF675616 is not Rana, it's a pathogen
		"Armillifer;" => 1, # tax issues
		"Meta;" => 1, # tax issues
		"Nippononeta;" => 1, # tax issues
		"Amanita;" => 1, # tax issues
		"Gymnophrys;" => 1, # tax issues
		"Malkaridae;" => 1, # tax issues
		"Tekelloides;" => 1, # tax issues
		"Tethya;" => 1, # tax issues
		"Aleiodes;" => 1, # tax issues
		"Chrysaora;" => 1, # tax issues
		"Choanocystis;" => 1,
		"Haploporus;" => 1,
		"Tritonus;" => 1, # KU711851 is a bettle not an arachnid
		"Physalacria;" => 1, # tax issues
		"Dracunculus;" => 1, 
		"Calyx;" => 1, # tax issues
		"Stramenopile;" => 1, # tax issues
		"Encalypta;" => 1, # tax issues
		"Tricladium;" => 1, # tax issues
		"Aquavolonida;" => 1, # tax issues
		"Halichondria;" => 1, # tax issues
		"Calliacantha;" => 1, # tax issues
		"Pterocystis;" => 1, # tax issues
		"Septobasidium;" => 1, # tax issues
		"Liriope;" => 1, 
		"Bumilleriopsis;" => 1,
		"Chaeteessa;" => 1,
		"Kerria;" => 1,
		"Paraglomus;" => 1, # tax issues
		"Hartaetosiga;" => 1, # tax issues
		"Hyalosira;" => 1, # tax issues
		"Hartaetosiga;" => 1, # tax issues
		"Haliplectus;" => 1, # tax issues
		"Oikomonas;" => 1, # tax issues
		"Dufourea;" => 1, # tax issues
		"Xanthonema;" => 1, # tax issues
		"Chaetophora;" => 1, 
		"Moissonia;" => 1, # GU194642 is a bug not an arachnid, tax issues
		"Florenciella;" => 1, # tax issues
		"Sternaspis;" => 1, # tax issues
		"Myxosporea;" => 1, # tax issues
		"Rissoella;" => 1, 
		"Pottia;" => 1, # tax issues
		"Rhinocladiella;" => 1, # tax issues
		"Tritonia;" => 1, # AJ224785 is oddball, tax issues here
		"Nephasoma;" => 1, # JN865022 is a worm not an arachnid
		"Serpula;" => 1,
		"Sarcinomyces;" => 1, # tax issues
		"Hydrachna;" => 1,
		"Moniliformis;" => 1, # tax issues
		"Dunaliella;" => 1,
#		"Pinus;" => 1,
		"Tubulicrinis;" => 1, # tax issues
		"Protaspis;" => 1, # tax issues
		"Ventrifissura;" => 1, # tax issues
		"Microcotylidae;" => 1, # tax issues
		"Palmophyllum;" => 1, # tax issues
		"Cirrophorus;" => 1, # tax issues
		"Corynaea;" => 1,
		"Castanea;" => 1,
		"Capnobotryella;" => 1, # tax issues
		"Massisteria;" => 1, # tax issues
		"Fonsecaea;" => 1,
		"Apochrysa;" => 1, # tax issues
		"Zantedeschia;" => 1,
		"Buxus;" => 1,
		"Myrothecium;" => 1, # taxonomic problems here
		"Glomus;" => 1, # tax issues
		"Chlorophyta;" => 1, # tax issues
		"Gnetum;" => 1,
		"Peltospira;" => 1, # AY923893 is not an arthropod
		"Grateloupia;" => 1,
		"Pentatomidae;" => 1, # GU194654 is a bug not an arachnid
		"Tetraspora;" => 1,
		"Thalassiosira;" => 1,
		"Bos;" => 1, # AAFC05032335 and AAFC05032349 are cattle, not fish
		"Hydnora;" => 1,
		"Phacotus;" => 1,
		"Melanthalia;" => 1,
		"Cytinus;" => 1,
		"Cymbomonas;" => 1,
		"Naineris;" => 1,
		"Sphaerospora;" => 1,
		"Posidonia;" => 1,
		"Gelidium;" => 1,
		"Facetotecta;" => 1,
		"Glycera;" => 1, # KT989347, KT989350, KT989341, and KT989346 are annelids not arthropods; GART01001051 looks like contaminant
		"Amphora;" => 1,
		"Phoenix;" => 1,
		"Diplogaster;" => 1, # tax issues
		"Arachis;" => 1, # tax issues
		"Mastocarpus;" => 1,
		"Pylaiella;" => 1,
		"Dollfusiella;" => 1, # tax issues
		"Monomorphina;" => 1,
		"Scenedesmus;" => 1,
		"Letharia;" => 1,
		"Orthocephalus;" => 1, # GU194647 is a bug not an arachnid
		"Cylindrotheca;" => 1,
		"Pandora;" => 1,
		"Cyanophora;" => 1,
		"Geodia;" => 1, # tax issues
		"Truncatella;" => 1,
		"Labrys;" => 1,
		"Europiella;" => 1, # GU194632 is a bug not an arachnid
		"Digenea;" => 1,
		"Torrubiella;" => 1, # tax issues
		"Humidicutis;" => 1, # tax issues
		"Polymerus;" => 1, # GU194667 is a bug not an arachnid
		"Oblongichytrium;" => 1, # tax issues
		"Bryopsis;" => 1,
		"Eucalyptus;" => 1, 
		"Dactylopodola;" => 1, # tax issues
		"Peltula;" => 1, # tax issues
		"Stemonitis;" => 1, # AY187085 has tax issue and JN123463 is an insect not a protist
		"Stenotus;" => 1, # GU194680 is a bug not an arachnid
		"Myiomma;" => 1, # GU194645 is bug not an arachnid
		"Adelphocoris;" => 1, # GU194607 is a bug not an arachnid
		"Cladonia;" => 1, # tax issues
		"Iris;" => 1,
		"Monalocoris;" => 1, # GU194643 is a bug not an arachnid
		"Prasinococcus;" => 1,
		"Junceella;" => 1, # AY962535 and AY962533 are not cnidarians
		"Chloropicon;" => 1,
		"Cimex;" => 1, # JQ782781 and JQ782778 are bugs not arachnids
		"Koshicola;" => 1,
		"Lichinella;" => 1, # tax issues
		"Anaplecta;" => 1, # KF855829 is oddball, tax issue
		"Latindia;" => 1, # tax issue with KP986332
		"Metarhizium;" => 1, # tax issues
		"Ostreobium;" => 1,
		"Albatrellus;" => 1, # tax issues
		"Nannochloropsis;" => 1, # tax issues between stramenopiles and chloroplastids
		"Metacrinus;" => 1, # KC626752 is an echinoderm not an arachnid
		"Macrolophus;" => 1, # GU194641 is a bug not an arachnid
		"Loxosomella;" => 1, # GU125746 is entoprocta not arthropoda
		"Prodesmodora;" => 1, # tax issues
		"Mertensiidae;" => 1, # MF599308 and MF599306 may be contaminants, tax issues here
		"Anomoloma;" => 1, # tax issues
		"Paracoccidioides;" => 1,
		"Bulboplastis;" => 1, 
		"Reticulamoeba;" => 1, # tax issues
		"Pertusaria;" => 1, # tax issues
		"Powellomyces;" => 1, # tax issues
		"Tuxedo;" => 1, # AY252328 is a bug not an arachnid
		"Pilophorus;" => 1,
		"Pigoraptor;" => 1, # tax issues
		"Beroe;" => 1, # tax issues
		"Pseudopediastrum;" => 1,
		"Diplomitoporus;" => 1, # tax issues
		"Orthotylus;" => 1, # GU194649, GU194652, GU194650 are bugs not arachnids
		"Olavius;" => 1, # AM233526 and AJ620508 are endosymbionts
		"Deraeocoris;" => 1, # tax issues
		"Lecanora;" => 1,
		"Sciaphila;" => 1,
		"Poria;" => 1,
		"Mitrastemon;" => 1,
		"Alaria;" => 1,
		"Stephanopyxis;" => 1, # tax issues
		"Bryocamptus;" => 1, # tax issues
		"Haliotrema;" => 1, # tax issues
		"Renouxia;" => 1,
		"Ganoderma;" => 1, # tax issues
		"Umbilicaria;" => 1, # tax issues
		"Mucor;" => 1, # tax issues
		"Veluticeps;" => 1, # tax issues
		"Tapinella;" => 1,
		"Rhodogorgon;" => 1,
		"Pleurotus;" => 1, # tax issues
		"Taiwania;" => 1, 
		"Italochrysa;" => 1, # tax issues
		"Sialis;" => 1, # EU815285 is a nematode contaminant, not an insect
		"Termitomyces;" => 1, # tax issues
		"Tricholomataceae;" => 1, # tax issues
		"Gigasporaceae;" => 1, # tax issues
		"Bothriocephalus;" => 1, # tax issues
		"Characiochloris;" => 1,
		"Acoela;" => 1, # tax issues
		"Trichaptum;" => 1, # tax issues
		"Echinorhynchida;" => 1, # tax issues
		"Oxyporus;" => 1,
		"Monomastix;" => 1,
		"Entomophthora;" => 1, # GEND01024600, GEND01024602, GENC01023619 are insect contaminants
		"Trigonostomum;" => 1,
		"Reclinomonas;" => 1,
		"Capsus;" => 1, # GU194620 is a bug, not an arachnid
		"Tubulipora;" => 1, # EU650325 is probably a contaminant
		"Jakoba;" => 1,
		"Heterosigma;" => 1,
		"Gadilida;" => 1, # tax issues
		"Incisitermes;" => 1, # FM160646 and AB032218 are symbionts NOT termites
		"Malawimonas;" => 1,
		"Corallina;" => 1,
		"Australomimetus;" => 1, # tax issues
		"Coelastrella;" => 1, # tax issues
		"Phaseolus;" => 1,
		"Mimetus;" => 1, # tax issues
		"Ero;" => 1, # tax issues
		"Polystilifera;" => 1, #tax issues
		"Lepidostroma;" => 1, # tax issues
		"Dioscorea;" => 1,
		"Kionochaeta;" => 1, # taxonomic issues
		"Kappaphycus;" => 1,
		"Marphysa;" => 1, # taxonomic issues here # AY040695 is an oddball
		"Epipyxis;" => 1,
		"Synarthrophyton;" => 1,
		"Aegisthidae;" => 1, # taxonomic issues here
		"Slimacomyces;" => 1, # taxonomic issues here
		"Cladochytrium;" => 1, # taxonomic issues here
		"Colossendeis;" => 1, # taxonomic issues here
		"Edwardsiella;" => 1,
		"Carteria;" => 1,
		"Euglena;" => 1,
		"Acrobeloides;" => 1, # taxonomic issues here
		"Entransia;" => 1,
		"Xylocoris;" => 1, # JQ782790 is an oddball, not an arachnid, it's a bug
		"Amastigomonas;" => 1, # taxonomic issues here
		"Apusomonadidae;" => 1,
		"Galeola;" => 1, 
		"Micractinium;" => 1,
		"Cassytha;" => 1,
		"Sebdenia;" => 1,
		"Frontonia;" => 1, # taxonomic issues, # FJ876952 is an oddball
		"Marsupiomonas;" => 1,
		"Strigamia;" => 1, # AFFK01003480 is a bacterial contaminant, not an arthropod
		"Euastacus;" => 1, # FJ965996 is probably a crustacean not an arachnid
		"Discosoma;" => 1, # taxonomic issues, # LT631296 is oddball
		"Aulactinia;" => 1, # taxonomic issues, # KT852129 is oddball
		"Nemertean;" => 1, # use of both Anopla and Enopla in SILVA taxonomy
		"Microglena;" => 1,
		"Phialemonium;" => 1, # taxonomic issues here, AB278185 is an oddball
		"Drosophila;" => 1, # EU188735 is an oddball, contaminant? taxonomic issue
		"Lotus;" => 1,
		"Asterionella;" => 1, # taxonomic issues here
		"Auxenochlorella;" => 1,
		"Lecythium;" => 1, # taxonomic issues here
		"Anodontia;" => 1, # taxonomic issues here
		"Andalucia;" => 1,
		"Gossypium;" => 1,
		"Fusarium;" => 1,
		"Dumontia;" => 1,
		"Vaucheria;" => 1,
		"Cryptobia;" => 1, # taxonomic issues
		"Gibellula;" => 1, # taxonomic issues
		"Fragaria;" => 1,
		"Paecilomyces;" => 1, # taxonomic issues
		"Vitis;" => 1,
		"Sarocladium;" => 1,
		"Chlorosarcina;" => 1,
		"Hirsutella;" => 1,
		"Spirogyra;" => 1,
		"Ophiocordyceps;" => 1, # taxonomic issues
		"Salpingoeca;" => 1, # MH490950 is an oddball, taxonomic issue
		"Cyanidium;" => 1,
		"Peridinium;" => 1, # taxonomic issues
		"Parachlorella;" => 1,
		"Nannochloris;" => 1, # taxonomic issues
		"Paramoebidium;" => 1, # taxonomic issues here
		"Hafniomonas;" => 1,
		"Glaucosphaera;" => 1,
		"Pycnococcus;" => 1,
		"Palpitomonas;" => 1,
		"Myristica;" => 1,
		"Nephroselmis;" => 1, 
		"Rhodella;" => 1,
		"Cordyceps;" => 1,
		"Colletotrichum;" => 1,
		"Mychonastes;" => 1, # AF106074 is an oddball, taxonomic issue
		"Poterioochromonas;" => 1,
		"Pterosperma;" => 1,
		"Fragilariaceae;" => 1, # JX413542 & JX413543 have taxonomic issues
		"Colacium;" => 1,
		"Blastodinium;" => 1, # JN257681 is an oddball, taxonomic issue
		"Prasinoderma;" => 1,
		"Hyalosynedra;" => 1, # MG684367 is an oddball, taxonomic issue
		"Lupinus;" => 1,
		"Grania;" => 1,
		"Tetradesmus;" => 1, # taxonomic problem # KX495035 is an oddball
		"Drosera;" => 1,
		"Adenoides;" => 1, # KX000292 is oddball, taxonomic issue
		"Botrydiopsis;" => 1, # taxonomic issues, # AJ579337 is an oddball
		"Cryptococcus;" => 1, # taxonomic issues with Filobasidiaceae/Cryptococcaceae
		"Tolypocladium;" => 1, # taxonomic issues here
		"Panax;" => 1, # taxonomic issues here GDQW01072153, GDQW01072150, and GDQW01072152 are oddballs
		"Neodangemannia;" => 1,
		"Selenidium;" => 1, # MF882905 & AY196708 are oddballs, taxonomic issues here
		"Closterium;" => 1, 
		"Dicyemennea;" => 1, # taxonomic issues
		"Treubia;" => 1,
		"Aecidium;" => 1, # taxonomic issues
		"Pinctada;" => 1, # taxonomic issues here
		"Spumella;" => 1, # taxonomic issues here
		"Corynactis;" => 1, # taxonomic issues here
		"Achnanthes;" => 1, # AJ535151 is an oddball, taxonomic problem here
		"Pyramimonas;" => 1,
		"Cymbidium;" => 1,
		"Strombomonas;" => 1,
		"Nakazawaea;" => 1, # taxonomic issues here
		"Oogamochlamys;" => 1,
		"Graphium;" => 1, # taxonomic issues here
		"Skeletonema;" => 1, # KY322507 is an oddball, taxonomic problem here
		"Kirchneriella;" => 1,
		"Penicillium;" => 1, # taxonomic problems with Talaromyces/Penicillium genera names both used in SILVA taxonomy
		"Puccinia;" => 1, # AORT01040852 is probably an arachnid contaminant not a fungus
		"Arthroderma;" => 1,
		"Zosterodasys;" => 1, # KC832953 is an oddball, taxonomic problem here
		"Parengyodontium;" => 1,
		"Karamea;" => 1, # KT302221 is probably an arachnid, not a bug
		"Lyctocoris;" => 1, # JQ782786 is probably a bug, not an arachnid
		"Leptographium;" => 1,
		"Pseudogymnoascus;" => 1,
		"Clonostachys;" => 1,
		"Trachelomonas;" => 1, # KT304871 is odd, taxonomic problem here
		"Lepocinclis;" => 1,
		"Cryptoglena;" => 1,
		"Corallimorphus;" => 1, # KJ483027 is odd, taxonomic problem here
		"Phacus;" => 1,
		"Prototheca;" => 1,
		"Andreaea;" => 1,
		"Orius;" => 1, # JQ782789 is probably a bug not an arachnid
		"Araucaria;" => 1,
		"Ophiostoma;" => 1,
		"Herpotrichiellaceae;" => 1, # taxonomic issues with # EU090194, consider removing
		"Grosmannia;" => 1,
		"Talaromyces;" => 1,
		"Liriodendron;" => 1,
		"Aspergillus;" => 1,
		"Cucumis;" => 1,
		"Cryptomonas;" => 1,
		"Chara;" => 1,
		"Hildenbrandia;" => 1,
		"Koliella;" => 1,
		"Dolichomastix;" => 1,
		"Gordionus;" => 1, # taxonomic issue here
		"Anthocoris;" => 1, # JQ782774 should probably be a bug, not an arachnid
		"Epispathidium;" => 1, # taxonomic issues here
		"Lobochlamys;" => 1,
		"Conchoecia;" => 1, # AF363296 is probably a contaminant
		"Crustomastix;" => 1, # KY980137 oddball
		"Paradoxia;" => 1,
		"Bracteacoccus;" => 1, # KF144165 & JX127179 oddballs
		"Stratiotes;" => 1,
		"Gymnodinium;" => 1, # taxonomic issues
		"Caulerpa;" => 1, # AF479702 & KY387611 oddballs, taxonomic issues here
		"Chloromonas;" => 1, # AF514403 oddball, taxonomic issue
		"Gaultheria;" => 1,
		"Darwinella;" => 1,
		"Chlorokybus;" => 1,
		"Kumanoa;" => 1, # KC511080 is probably a contaminant/taxonomic issue here
		"Gracilaria;" => 1,
		"Lecanicephalidea;" => 1, # KF685783 & KF685792 oddballs, taxonomic issues here
		"Apicoporus;" => 1, # KP790148 oddball, taxonomic issue here
		"Chlorochytrium;" => 1, # HE687020 & HE860259 oddballs, taxonomic issues here
		"Asclepias;" => 1,
		"Proteocephalus;" => 1, # X99976 oddball, taxonomic issue here?
		"Ceramium;" => 1,
		"Panorpa;" => 1, # DQ008172 oddball, misidentification/taxonomic issue here?
		"Karlodinium;" => 1,
		"Siphonostomatoida;" => 1, # MF077749 oddball, taxonomic issues here?
		"Brassica;" => 1,
		"Magelona;" => 1, # U50969 oddball, taxonomic issues here?
		"Picochlorum;" => 1,
		"Vaccinium;" => 1,
		"Sporolithon;" => 1,
		"Hypsibius;" => 1,
		"Epiactis;" => 1, # Z92904 oddball, taxonomic issues here?
		"Allium;" => 1,
		"Acorus;" => 1,
		"Wisteria;" => 1,
		"Pavlova;" => 1, # DQ075199 is probably haptophyte, not a plant
		"Pseudochlorella;" => 1,
		"Oltmannsiellopsis;" => 1,
		"Alcyonidium;" => 1, # FJ196130 is probably a bryozoan, not arthropod
		"Virgularia;" => 1, # taxonomic issues here
		"Radix;" => 1, # MUZC01002578 oddball, contaminant?
		"Chlorodesmis;" => 1,
		"Coryne;" => 1, # GQ424325 oddball, taxonomic issues here?
		"Oryza;" => 1, 
		"Solifugae;" => 1, # KY573618 is an oddball, taxonomic issues here?
		"Halimeda;" => 1,
		"Herdmania;" => 1,
		"Plectus;" => 1, # MG993558 should probably be a nemotode too, not an arthropod
		"Neogoniolithon;" => 1,
		"Chaetosphaeridium;" => 1,
		"Botryococcus;" => 1, 
		"Ammonia;" => 1, # U07937 is oddball, taxonomic issues here
		"Parabothriocephalus;" => 1, # taxonomic issues here, not enough info
		"Philodina;" => 1, # AF154567 is an oddball, taxonomic issues?
		"Dicyema;" => 1, # LT669877 is oddball
		"Coccomyxa;" => 1,
		"Nerita;" => 1, # AY923889 may be a contaminant, its probably a gastropod not an arachnid
		"Thorea;" => 1,
		"Vigna;" => 1,
		"Chaetopterus;" => 1, # KX896469 is oddball, taxonomic issues here?
		"Pseudopleurococcus;" => 1, # KM020154 is an oddball, taxonomic issues here?
		"Phoronis;" => 1, # GU125758 is an oddball, taxonomic issues here?
		"Golenkinia;" => 1,
		"Hydra;" => 1,
		"Actinocyclus;" => 1,
		"Chlorogonium;" => 1,
		"Eimeria;" => 1,
		"Zygnema;" => 1,
		"Pedinomonas;" => 1,
		"Pedinophyceae;" => 1,
		"Plumatella;" => 1, # EU650324 is probably not an arthropod
		"Xylochloris;" => 1,
		"Berkeleya;" => 1, # uncertain taxonomy?
		"Plantago;" => 1,
		"Tribonema;" => 1, #problematic taxonomy?
		"Ulva;" => 1,
		"Amoebophrya;" => 1, # problematic taxonomy?
		"Thyasira;" => 1, # LC187038 is probably a fungal contaminant
		"Dictyocha;" => 1, # SAR inserted between Eukaryota and Stramenopiles sometimes
		"Kofoidinium;" => 1, #problematic taxonomy? not enough examples to know for sure
		"Pogonatum;" => 1, # AY126968 needs a look
		"Anthocephalum;" => 1, # KM658178 is probably a flatworm, not an arachnid
		"Acanthamoeba;" => 1,
		"Chloroidium;" => 1,
		"Chaetoceros;" => 1,
		"Paralia;" => 1, # AJ535174 needs a second look
		"Taenioma;" => 1,
		"Stichococcus;" => 1,
		"Fucus;" => 1,
		"Vannella;" => 1,
		"Jenufa;" => 1,
		"Arabidopsis;" => 1,
		"Amphidinium;" => 1, # GU295202 odd ball lineage
		"Chondria;" => 1,
		"Tetraselmis;" => 1,
		"Citrus;" => 1,
		"Acanthocephala;" => 1,
		"Ochromonas;" => 1, # EF165142 should probably be Ochromonas, not Uroglena; all the Ochromonadales should probably be Chromulinales
		"Oikopleura;" => 1,
		"Robinia;" => 1,
		"Ettlia;" => 1,
		"Solanum;" => 1,
		"Hordeum;" => 1,
		"Synura;" => 1,
		"Prorocentrum;" => 1,
		"Myrmecia;" => 1,
		"Pilostyles;" => 1,
		"Odontella;" => 1,
		"Sargassum;" => 1,
		"Klebsormidium;" => 1,
		"Thecadinium;" => 1,
		"Dracoderes;" => 1, # LC032113 should probably be Scalidophora like the rest
		"Spirotrichonympha;" => 1,
		"Heteronemertea;" => 1, # MF512067 and MF512069 should probably be Heteronemertea like the rest
		"Wildemania;" => 1,
		"Perideraion;" => 1, # MH063516 should probably be Fragilariales like the rest
		"Chondracanthus;" => 1,
		"Plagiogramma;" => 1,
		"Pygophora;" => 1,
		"Sycon;" => 1, # D15066 should probably be Leucosolenida like the rest
		"Coleochaete;" => 1,
		"Chlorella;" => 1,
		"Pediastrum;" => 1,
		"Linum;" => 1,
		"Azadinium;" => 1, # MF033111 should probably be Gonyaulacales like the rest
		"Haplosclerida;" => 1, # KX894474 should probably be Demospongiae like the rest
		"Ixodes;" => 1, # GBXQ01023944 and GBXQ01000923 are probably incorrectly annotated in SILVA, should be arachnids like the rest
		"Dactylanthus;" => 1,
		"Daphnia;" => 1, # Z23111 may be a parasite of Daphnia???
		"Bangia;" => 1,
		"Gynostemma;" => 1,
		"Cynomorium;" => 1,
		"Medicago;" => 1, 
		"Dinophyceae;" => 1, 
		"Pyropia;" => 1, 
		"Nitzschia;" => 1, 
		"Ostreococcus;" => 1,
		"Craniella;" => 1, # KY652840 is an oddball with different taxonomic lineage
		"Plasmodium;" => 1, 
		"Babesia;" => 1, 
		"Apicomplexa;" => 1, 
		"Labyrinthulomycetes;" => 1, # because of Amorphea;Amoebozoa between Eukaryota & SAR WTF.
		"Coptosoma;" => 1, # KJ461259.1.1742 may be misclassified (GenBank says Hemiptera, not Coleoptera)
		"Chlamydomonas;" => 1, 
		"Poecilosclerida;" => 1,
		"Zea;" => 1, 
		"Porphyra;" => 1,
		"basidiomycete;" => 1,
		"marine_group;" => 1,
		"undef;" => 1,
		"unidentified;" => 1,
		"uncultured;" => 1,
		"Incertae_Sedis;" => 1,
		"Candida;" => 1,
		"Rhizophagus;" => 1,
		"Navicula;" => 1,
		"Karenia;" => 1,
		"Nectria;" => 1,
		"Lobophora;" => 1,
		"Zonaria;" => 1,
		"Ctenophora;" => 1,
		"Campanella;" => 1,
		"Uronema;" => 1,
		"Lessonia;" => 1,
		"Atractomorpha;" => 1,
		"Bostrychia;" => 1,
		"Pontogeneia;" => 1,
		"Mastophora;" => 1,
		"Protura;" => 1,
		"Parachela;" => 1,
		"Acrotylus;" => 1,
		"Plecoptera;" => 1,
		"Achlya;" => 1,
		"Nemastoma;" => 1,
		"Choanoflagellida;" => 1,
		"Ptilophora;" => 1,
		"Olea;" => 1,
		"Pterospora;" => 1,
		"Dacrydium;" => 1,
		"Beauveria;" => 1,
		"Phoma;" => 1,
		"Drymonia;" => 1,
		"Roya;" => 1,
		"Heterococcus;" => 1,
		"Actinobacteria;" => 1,
		"Rhodococcus;" => 1,
		"Bacillus;" => 1,
		"Family_XI;" => 1,
		"Unknown_Family;" => 1,
		"Thermotogae;" => 1,
		"Aquificae;" => 1,
		"Family_XII;" => 1,
		"Deferribacteres;" => 1,
		"Chlamydiae;" => 1,
		"Nitrospira;" => 1,
		"Latescibacteria;" => 1,
		"Chrysiogenetes;" => 1,
		"Gemmatimonadetes;" => 1,
		"Planococcus;" => 1,
		"Elusimicrobia;" => 1,
		"Paracoccus;" => 1);

#read in the SILVA taxonomy file (domain -> genus)
open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@tax = <IN>;
close IN;

#hash the taxonomy file 
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;
 
	@line = split(/\t/,$line);
	$path = $line[0];
	$rank = $line[2];

	$path =~ s/\s+/_/g;
	$path =~ s/\(//g;
	$path =~ s/\)//g;
	$path =~ s/\.//g;
	$path =~ s/\</_/g;
	$path =~ s/\>/_/g;

	$map{$path} = $rank;

	$i++;
}
$i=0;

#read in problem_species file that contains list of species with more than one taxonomic path
open (IN2, "<", $ARGV[2]) || die "Cannot open problem species infile: $!\n";
@problem = <IN2>;
close IN2;

#hash the problem species
while ($problem[$i]) {
	$line = $problem[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$species = $line[0];
	$lineage = $line[1];
	$abundance = $line[2];

	$species =~ s/\s+/_/g;
	$species =~ s/\(//g;
	$species =~ s/\)//g;
	$species =~ s/\.//g;
	$species =~ s/\</_/g;
	$species =~ s/\>/_/g;

	$lineage =~ s/\s+/_/g;
	$lineage =~ s/\(//g;
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g;
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;

	$problem{$species}{$lineage} = $abundance;
	$i++;
}
$i=0;

#print Dumper(\%problem);

# remove species from %problem if the dominant lineage can be salvaged into new %fixproblem
# loop through hash of hashes
# find inner key with largest value, indicate if multiple keys tied for largest value
foreach $species (keys %problem) {

	($newlineage, $prob) = key_with_highest_val($problem{$species});

	if ($prob == 0) { # largest value found, else prob == 1 indicates a tie
		$fixproblem{$species} = $newlineage; 
		delete($problem{$species}{$newlineage}); # don't delete whole species, just the most abundant lineage from %problem
	}
}

# read in the SILVA fasta.gz file
@fasta = `zcat $ARGV[1]`;

# create an outfile
open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";


# populate %duplicates
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	# account for duplicate genera that appeared when species names were labelled as 'Mitochondria' or 'Chloroplast' 
	if (($line =~ /Mitochondria;/) || ($line =~ /Chloroplast;/)) {
		print "found $line with mit or chlor\n"; # test
		my @splitline = split(';',$line);
		$species = pop(@splitline);
		my @species = split(' ',$species);
		my $genusFromSpecies = shift(@species);
		$genusFromSpecies = $genusFromSpecies.";";
		$duplicates{$genusFromSpecies}= 1;
		@splitline=();
		@species=();
		$genusFromSpecies=();
#		print Dumper(\%duplicates);
	}
	$i++;
}
$i=0;

print Dumper(\%duplicates);


# parse through SILVA fasta file
# check for higher level duplicate taxa that needs to be reformatted
# check for problem and fixproblem taxa
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	#parse FASTA header
	if ($line =~ /^>/) {
		if (length $seq) { # check if a sequence exists
			if (($euk == 1 && $prob == 0) || ($euk == 1 && $prob == 2)) {
				print OUT "$seq\n";
				$seq="";
			}
		}

		if ($line =~ /Eukaryota;/) { #just process Eukaryota for now
			$euk=1;

			parse_header();
			reformat_lineage();

			if ($prob == 0) {
				print OUT ">$acc cellularOrganisms;$newlineage\n";
				$good_seq_counter++;
			}
			elsif ($prob == 2) {
				parse_header2();
				reformat_lineage2();

				if ($prob == 1) {
					$i++;
					%observed=();
					$newlineage=();
					next;
				}

				print OUT ">$acc cellularOrganisms;$newlineage\n";
				$salvaged_seq_counter++;
			}
			else { # prob == 1 indicates species has multiple lineages in SILVA taxonomy file, can't choose a best one (tie among choices)
				$bad_seq_counter++;
				$i++;
				%observed=();
				$newlineage=();
				next;
			}
		}
		else {
			$euk=0;
			%observed=();
			$newlineage=();
			$i++;
			next;
		}
	}
	#concatenate sequence across lines
	elsif (($euk == 1 && $prob == 0) || ($euk == 1 && $prob == 2)) {
		$line =~ s/U/T/g; #convert U's to T's
		$line = lc $line; #convert to lower case
		$seq = $seq.$line; #append new line to old line
	}
	$i++;
	%observed=();
	$newlineage=();
}
$i=0;

#don't forget to print last seq if appropriate
if (($euk == 1 && $prob == 0) || ($euk == 1 && $prob == 2)) {
	print OUT "$seq\n";
}
$euk=();
$prob=();
$seq="";

# print status of seqs (good, salvaged, problems)
print "Eukrayota good: $good_seq_counter\n";
print "Eukaryota salvaged: $salvaged_seq_counter\n";
print "Eukaryota bad: $bad_seq_counter\n";

$good_seq_counter=0;
$salvaged_seq_counter=0;
$bad_seq_counter=0;

# read in the prokaryote silva.centroids.gz file
@centroids = `zcat $ARGV[3]`;

# parse through outgroup centroids fasta file
# check for higher level duplicate taxa that needs to be reformatted
# check for problem and fixproblem taxa
while ($centroids[$i]) {
	$line = $centroids[$i];
	chomp $line;

	#parse FASTA header
	if ($line =~ /^>/) {
		if (length $seq) { # check if a sequence exists
			if (($euk == 0 && $prob == 0) || ($euk == 0 && $prob == 2)) {
				print OUT "$seq\n";
				$seq="";
			}
		}

		if ($line =~ /Chloroplast;/) { # account for SSU sequence from a chloroplast
			$line = $line."_Chloroplast";
		}
		if ($line =~ /Mitochondria;/) { # account for SSU sequence from mitochondria
			$line = $line."_Mitochondrion";
		}

		if (($line =~ /Bacteria;/) || ($line =~ /Archaea;/)) { #just process Eukaryota for now
			$euk=0;

			parse_header_pro();
			reformat_lineage();

			if ($prob == 0) {
				print OUT ">$acc cellularOrganisms;$newlineage\n";
				$good_seq_counter++;
			}
			elsif ($prob == 2) {
				parse_header2_pro();
				reformat_lineage2();

				if ($prob == 1) {
					$i++;
					%observed=();
					$newlineage=();
					next;
				}

				print OUT ">$acc cellularOrganisms;$newlineage\n";
				$salvaged_seq_counter++;
			}
			else { # prob == 1 indicates species has multiple lineages in SILVA taxonomy file, can't choose a best one (tie among choices)
				$i++;
				%observed=();
				$newlineage=();
				$bad_seq_counter++;
				next;
			}
		}
		else {
			$euk=1;
			$i++;
			%observed=();
			$newlineage=();
			next;
		}
	}
	#concatenate sequence across lines
	elsif (($euk == 0 && $prob == 0) || ($euk == 0 && $prob == 2)) {
		$line =~ s/U/T/g; #convert U's to T's
		$line = lc $line; #convert to lower case
		$seq = $seq.$line; #append new line to old line
	}
	$i++;
	%observed=();
	$newlineage=();
}
$i=0;

#don't forget to print last seq if appropriate
if (($euk == 0 && $prob == 0) || ($euk == 0 && $prob == 2)) {
	print OUT "$seq\n";
}
close OUT;

# print status of seqs (good, salvaged, problems)
print "Prokaryote good: $good_seq_counter\n";
print "Prokaryote salvaged: $salvaged_seq_counter\n";
print "Prokaryote bad: $bad_seq_counter\n";

# run some commands to reshape the output
print `grep ">" testNBC.fasta | awk 'BEGIN { FS= " " }; { print \$2}'  >> taxa.awk`;
print `sed 's/;/\t/g' taxa.awk >> taxa.awk.sed`;
print `sort -u taxa.awk.sed >> taxa.awk.sed.uniq`;


##########

sub parse_header {
	
	%observed=();
	$line =~ s/^>//g;
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;

	@lineage = split(/;/,$lineage);
	$species = pop(@lineage);
	$observed{"species"} = $species;

	# reconstruct domain->genus lineage to effectively check %problem
	# no need for termianl ";"
	# keep unedited originalLineage to search in %problem 
	$originalLineage = join(';', @lineage);
#	print "species: $species originalLineage: $originalLineage\n";

	if (($species =~ /^[A-Z]/) && ($species !~ /^Eukaryota/) 
			&& ($species !~ /^Arthropoda/) 
			&& ($species !~ /^Tylenchida/)
			&& ($species !~ /^Cercozoa/)
			&& ($species !~ /^Polychaeta/)
			&& ($species !~ /^Ptychoderidae/)
			&& ($species !~ /^Palaeonemertea/)
			&& ($species !~ /^Trebouxiophyceae/)
			&& ($species !~ /^Chlorellales/)
			&& ($species !~ /^Chlorellaceae/)
			&& ($species !~ /^Peridiniales/)
			&& ($species !~ /^Harrimaniidae/)
			&& ($species !~ /^Ulvophyceae/)
			&& ($species !~ /^Dinobryon/)
			&& ($species !~ /^Heteronemertea/)
			&& ($species !~ /^Reticulitermes/)
			&& ($species !~ /^Cryptotermes/)
			&& ($species !~ /^Neotermes/)
			&& ($species !~ /^Myxozoa/)
			&& ($species !~ /^Solemya/)
			&& ($species !~ /^endosymbiont/)
			&& ($species !~ /^Rhinebothriidea/)
			&& ($species !~ /^Jakobida/)
			&& ($species !~ /^Nematoda/)
			&& ($species !~ /^Xanthophyceae/)
			&& ($species !~ /^Monostilifera/)
			&& ($species !~ /^Scenedesmaceae/)
			&& ($species !~ /^Tubulinea/)
			&& ($species !~ /^Sordariomycetes/)
			&& ($species !~ /^Foraminifera/)
			&& ($species !~ /^Psocoptera/)
			&& ($species !~ /^Siphonostomatoida/)
			&& ($species !~ /^Tetraphyllidea/)
			&& ($species !~ /^Eurotiomycetes/)
			&& ($species !~ /^Chlorophyceae/)
			&& ($species !~ /^Microascaceae/)
			&& ($species !~ /^Haptoria/)
			&& ($species !~ /^Pleosporales/)
			&& ($species !~ /^Phaeothamniophyceae/)
			&& ($species !~ /^Enteropneusta/)
			&& ($species !~ /^Tubulanidae/)
			&& ($species !~ /^Diptera/)
			&& ($species !~ /^Tricladida/)
			&& ($species !~ /^Amoebozoa/)
			&& ($species !~ /^Bicosoecida/)
			&& ($species !~ /^Ectocarpales/)
			&& ($species !~ /^Pauropoda/)
			&& ($species !~ /^Pelagophyceae/)
		 	&& ($species !~ /^Chrysophyceae/)
			&& ($species !~ /^Eustigmatophyceae/)
			&& ($species !~ /^Thraustochytriidae/)
			&& ($species !~ /^Dactylogyridae/)
			&& ($species !~ /^Lasiosphaeriaceae/)
			&& ($species !~ /^Tardigrada/)
			&& ($species !~ /^Bdelloidea/)
			&& ($species !~ /^Crinoidea/)
			&& ($species !~ /^Rhizaria/)
			&& ($species !~ /^Spumella-like/)
			&& ($species !~ /^Siphonophorae/)
			&& ($species !~ /^Coccidia/)
			&& ($species !~ /^Chytridiales/)
			&& ($species !~ /^Roussoellaceae/)) { # A real genus name should be capitalized, else skip this part
	
		my @species = split('_',$species); # do this to accomodate for genera missing from SILVA tax file ex. Zea
		my $genusFromSpecies = shift(@species);
		my $genusFromLineage = pop(@lineage);

		if ($genusFromSpecies eq $genusFromLineage) {
			$observed{"genus"} = $genusFromLineage;
		}
		else {
			$observed{"genus"} = $genusFromSpecies;
			my $scalar = scalar(@lineage);
			$lineage[$scalar] = $genusFromLineage;
		}
	}

	foreach $taxon (@lineage) {
		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
		@taxon = split(/;/,$fulltaxon);
		$taxon = pop(@taxon); # use this to add to %observed

		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
			$rank = $map{$fulltaxon};
			$observed{$rank} = $taxon;
		}
		$prevtaxon = $fulltaxon;
	}
	$prevtaxon="";
}

##########

sub reformat_lineage {
	if (exists $observed{"domain"}) {
		$domain = $observed{"domain"};
		$domain = $domain.";";
	}
	else {
		print "Couldn't find domain for accession $acc from reformat_lineage\n";
		$domain="unknownDomain;";
	}

	if (exists $observed{"kingdom"}) {
		$kingdom = $observed{"kingdom"};
		$kingdom = $kingdom.";";
	
		if (exists $duplicates{$kingdom}) {
			$domain =~ s/;$//;
			$kingdom = $domain."_".$kingdom;
			$domain = $domain.";";
		}
	}
	else {
		$domain =~ s/;$//g;
		$kingdom = $domain."_undef;";
		$domain = $domain.";";
	}

	if (exists $observed{"phylum"}) {
		$phylum = $observed{"phylum"};
		$phylum = $phylum.";";

		if (exists $duplicates{$phylum}) {# handle duplicate phyla
			$kingdom =~ s/;$//;
			$phylum = $kingdom."_".$phylum;
			$kingdom = $kingdom.";";
		}
	}
	else {
		$kingdom =~ s/;$//;
		$phylum = $kingdom."_undef;";
		$kingdom = $kingdom.";";
	}

	if (exists $observed{"class"}) {
		$class = $observed{"class"};
		$class = $class.";";

		if (exists $duplicates{$class}) { #handle duplicate classes
			$phylum =~ s/;$//;
			$class = $phylum."_".$class;
			$phylum = $phylum.";";
		}
	}
	else {
		$phylum =~ s/;$//;
		$class = $phylum."_undef;";
		$phylum = $phylum.";";
	}

	if (exists $observed{"order"}) {
		$order = $observed{"order"};
		$order = $order.";";

		if (exists $duplicates{$order}) { #handle duplicate orders
			$class =~ s/;$//;
			$order = $class."_".$order;
			$class = $class.";";
		}
	}
	else {
		$class =~ s/;$//;
		$order = $class."_undef;";
		$class = $class.";";
	}

	if (exists $observed{"family"}) {
		$family = $observed{"family"};
		$family = $family.";";

		if (exists $duplicates{$family}) {#handle duplicate families
			$order =~ s/;$//;
			$family = $order."_".$family;
			$order = $order.";";
		}
	}
	else {
		$order =~ s/;$//;
		$family = $order."_undef;";
		$order = $order.";";
	}

	if (exists $observed{"genus"}) {
		$genus = $observed{"genus"};
		$genus = $genus.";";

		if (exists $duplicates{$genus}) {# handle duplicate genera
			$family =~ s/;$//;
			$genus = $family."_".$genus;
			$family = $family.";";
		}
	}
	else {
		$family =~ s/;$//;
		$genus = $family."_undef;";
		$family = $family.";";
#		$prob = 1;
#		goto end;
	}

	if (exists $observed{"species"}) {
		$species = $observed{"species"};

#		if (exists $duplicates{$species}) { # handle duplicate species, "unidentified"
#			print "Found $species in %duplicates\n";
#			$prob=1;
#			goto end;
#		}
		
#		print "species: ".$species."original lineage :".$originalLineage."\n";
		if (exists $problem{$species} && exists $problem{$species}{$originalLineage}) { #don't process species with more than one taxonomic path
			$prob=1;
			goto end;
		}
		
		elsif (exists $fixproblem{$species}) {
			my $templineage = $fixproblem{$species};
			$prob=2;
			goto end;
		}
		else {
			$prob=0;
			$newlineage = $domain.$kingdom.$phylum.$class.$order.$family.$genus.$species;
		}
	}
	else {
		$prob = 1;
		goto end;
	}

	$domain=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();

	end:
	$domain=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
	return $prob;

}

##########

sub parse_header2 {
#	print "parse_header2\n";
	%observed=();
	$line =~ s/^>//g;
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;

	# parse the original lineage
	@lineage = split(/;/,$lineage);
	$species = pop(@lineage);
	$observed{"species"} = $species;

	# most abundant lineage from problem species
	$lineage = $fixproblem{$species};
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;
#	print "lineage from parse_header2: $lineage\n";

	@lineage = split(/;/, $lineage); # domain -> genus

	# Do this to accomodate for genera missing from SILVA tax file ex. Zea
	if (($species =~ /^[A-Z]/) 
			&& ($species !~ /^Eukaryota/) # A real species would begin with a capital, else skip this part
			&& ($species !~ /^Arthropoda/) # These capitalized naems are not the genus part of a Latin binomial
			&& ($species !~ /^Tylenchida/)
			&& ($species !~ /^Cercozoa/)
			&& ($species !~ /^Polychaeta/)
			&& ($species !~ /^Ptychoderidae/)
			&& ($species !~ /^Palaeonemertea/)
			&& ($species !~ /^Trebouxiophyceae/)
			&& ($species !~ /^Chlorellales/)
			&& ($species !~ /^Chlorellaceae/)
			&& ($species !~ /^Peridiniales/)
			&& ($species !~ /^Harrimaniidae/)
			&& ($species !~ /^Ulvophyceae/)
			&& ($species !~ /^Dinobryon/)
			&& ($species !~ /^Heteronemertea/)
			&& ($species !~ /^Reticulitermes/)
			&& ($species !~ /^Cryptotermes/)
			&& ($species !~ /^Neotermes/)
			&& ($species !~ /^Myxozoa/)
			&& ($species !~ /^Solemya/)
			&& ($species !~ /^endosymbiont/)
			&& ($species !~ /^Rhinebothriidea/)
			&& ($species !~ /^Jakobida/)
			&& ($species !~ /^Nematoda/)
			&& ($species !~ /^Xanthophyceae/)
			&& ($species !~ /^Monostilifera/)
			&& ($species !~ /^Scenedesmaceae/)
			&& ($species !~ /^Tubulinea/)
			&& ($species !~ /^Sordariomycetes/)
			&& ($species !~ /^Foraminifera/)
			&& ($species !~ /^Psocoptera/)
			&& ($species !~ /^Siphonostomatoida/)
			&& ($species !~ /^Tetraphyllidea/)
			&& ($species !~ /^Eurotiomycetes/)
			&& ($species !~ /^Chlorophyceae/)
			&& ($species !~ /^Microascaceae/)
			&& ($species !~ /^Haptoria/)
			&& ($species !~ /^Pleosporales/)
			&& ($species !~ /^Phaeothamniophyceae/)
			&& ($species !~ /^Enteropneusta/)
			&& ($species !~ /^Tubulanidae/)
			&& ($species !~ /^Diptera/)
			&& ($species !~ /^Tricladida/)
			&& ($species !~ /^Amoebozoa/)
			&& ($species !~ /^Bicosoecida/)
			&& ($species !~ /^Ectocarpales/)
			&& ($species !~ /^Pauropoda/)
			&& ($species !~ /^Pelagophyceae/)
			&& ($species !~ /^Chrysophyceae/)
			&& ($species !~ /^Eustigmatophyceae/)
			&& ($species !~ /^Thraustochytriidae/)
			&& ($species !~ /^Dactylogyridae/)
			&& ($species !~ /^Lasiosphaeriaceae/)
			&& ($species !~ /^Tardigrada/)
			&& ($species !~ /^Bdelloidea/)
			&& ($species !~ /^Crinoidea/)
			&& ($species !~ /^Rhizaria/)
			&& ($species !~ /^Spumella-like/)
			&& ($species !~ /^Siphonophorae/)
			&& ($species !~ /^Coccidia/)
			&& ($species !~ /^Chytridiales/)
			&& ($species !~ /^Roussoellaceae/)) {
	
		my @species = split('_',$species); 
		my $genusFromSpecies = shift(@species);
		my $genusFromLineage = pop(@lineage);

		if ($genusFromSpecies eq $genusFromLineage) {
			$observed{"genus"} = $genusFromLineage;
		}
		else {
			$observed{"genus"} = $genusFromSpecies;
			my $scalar = scalar(@lineage);
			$lineage[$scalar] = $genusFromLineage;
		}
#		print "parse_header2 lineage @lineage \n"; #testing
	}

	# parse the most abundant lineage 
	foreach $taxon (@lineage) {
		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
#		print "fulltaxon: $fulltaxon\n";
		@taxon = split(/;/,$fulltaxon);
		$taxon = pop(@taxon); # use this to add to %observed

		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
			$rank = $map{$fulltaxon};
			$observed{$rank} = $taxon;
#			print "parse_header2: rank $rank taxon $taxon\n";
		}
		$prevtaxon = $fulltaxon;
	}
#	print "parse_header2 %observed\n";
#	print Dumper(\%observed);
	$prevtaxon="";
}

##########

sub reformat_lineage2 {
#	print "reformat_lineage2\n";
	if (exists $observed{"domain"}) {
		$domain = $observed{"domain"};
		$domain = $domain.";";
	}
	else {
		print "Couldn't find domain for accession $acc in reformat_lineage2\n";
		$domain="unknownDomain;";
	}

	if (exists $observed{"kingdom"}) {
		$kingdom = $observed{"kingdom"};
		$kingdom = $kingdom.";";

		if (exists $duplicates{$kingdom}) {
			$domain =~ s/;$//;
			$kingdom = $domain."_".$kingdom;
			$domain = $domain.";";
		}
	}
	else {
		$domain =~ s/;$//g;
		$kingdom = $domain."_undef;";
		$domain = $domain.";";
	}

	if (exists $observed{"phylum"}) {
		$phylum = $observed{"phylum"};
		$phylum = $phylum.";";

		if (exists $duplicates{$phylum}) {# handle duplicate phyla
			$kingdom =~ s/;$//;
			$phylum = $kingdom."_".$phylum;
			$kingdom = $kingdom.";";
		}
	}
	else {
		$kingdom =~ s/;$//;
		$phylum = $kingdom."_undef;";
		$kingdom = $kingdom.";";
	}

	if (exists $observed{"class"}) {
		$class = $observed{"class"};
		$class = $class.";";

		if (exists $duplicates{$class}) { #handle duplicate classes
			$phylum =~ s/;$//;
			$class = $phylum."_".$class;
			$phylum = $phylum.";";
		}
	}
	else {
		$phylum =~ s/;$//;
		$class = $phylum."_undef;";
		$phylum = $phylum.";";
	}

	if (exists $observed{"order"}) {
		$order = $observed{"order"};
		$order = $order.";";

		if (exists $duplicates{$order}) { #handle duplicate orders
			$class =~ s/;$//;
			$order = $class."_".$order;
			$class = $class.";";
		}
	}
	else {
		$class =~ s/;$//;
		$order = $class."_undef;";
		$class = $class.";";
	}

	if (exists $observed{"family"}) {
		$family = $observed{"family"};
		$family = $family.";";

		if (exists $duplicates{$family}) {#handle duplicate families
			$order =~ s/;$//;
			$family = $order."_".$family;
			$order = $order.";";
		}
	}
	else {
		$order =~ s/;$//;
		$family = $order."_undef;";
		$order = $order.";";
	}

	if (exists $observed{"genus"}) {
		$genus = $observed{"genus"};
		$genus = $genus.";";

		if (exists $duplicates{$genus}) { # handle duplicate genera
			$family =~ s/;$//;
			$genus = $family."_".$genus;
			$family = $family.";";
		}
	}
	else {
		$prob = 1;
		goto end;
	}

	if (exists $observed{"species"}) {
		$species = $observed{"species"};

		if (exists $duplicates{$species}) { # this shouldn't even happen
#			print "Found $species in %duplicates\n";
			$prob=1;
			goto end;
		}
		else {
			$newlineage = $domain.$kingdom.$phylum.$class.$order.$family.$genus.$species;
#			print "new lineage: $newlineage\n";
		}
	}
	else {
		$prob = 1;
#		print "couldn't find $species in %observed\n";
		goto end;
	}	

	$domain=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();

	end:
	$domain=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();

		return $prob;

}

##########
# find the inner key of a hash of hashes with the largest value
# keep track of most abundant lineage; track minor lineages; report ties
# https://stackoverflow.com/questions/35401306/find-key-for-greatest-value-in-hash-of-hashes-in-perl
sub key_with_highest_val {

	my ($h) = @_; #hash of hashes where outer key is species and already specified
	my $hi_v; # value with highest abundance
	my $hi_k; # innerkey with highest abundance aka lineage
	my $tie=0;

	# loop through each lineage, outer keys
	for my $k (keys(%$h)) {
		my $v = $h->{$k};
		if (!defined($hi_v) || $v > $hi_v) {
			$hi_v = $v;
			$hi_k = $k;
			$tie = 0;
		}
		elsif ($v == $hi_v) {
			$tie = 1;
		}
		else {
			$tie = 0;
		}
	}

	return ($hi_k, $tie);
}

##########

sub parse_header_pro {
	%observed=();
	$line =~ s/^>//g;
	$line =~ s/_/ /; # replace FIRST underscore with space
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;

	@lineage = split(/;/,$lineage);

	$species = pop(@lineage);
	$observed{"species"} = $species;
	
	if (($species =~ /^[A-Z]/) && ($species !~ /^Eukaryota/) 
			&& ($species !~ /^Arthropoda/) 
			&& ($species !~ /^Tylenchida/)
			&& ($species !~ /^Cercozoa/)
			&& ($species !~ /^Polychaeta/)
			&& ($species !~ /^Ptychoderidae/)
			&& ($species !~ /^Palaeonemertea/)
			&& ($species !~ /^Trebouxiophyceae/)
			&& ($species !~ /^Chlorellales/)
			&& ($species !~ /^Chlorellaceae/)
			&& ($species !~ /^Peridiniales/)
			&& ($species !~ /^Harrimaniidae/)
			&& ($species !~ /^Ulvophyceae/)
			&& ($species !~ /^Dinobryon/)
			&& ($species !~ /^Heteronemertea/)
			&& ($species !~ /^Reticulitermes/)
			&& ($species !~ /^Cryptotermes/)
			&& ($species !~ /^Neotermes/)
			&& ($species !~ /^Myxozoa/)
			&& ($species !~ /^Solemya/)
			&& ($species !~ /^endosymbiont/)
			&& ($species !~ /^Rhinebothriidea/)
			&& ($species !~ /^Jakobida/)
			&& ($species !~ /^Nematoda/)
			&& ($species !~ /^Xanthophyceae/)
			&& ($species !~ /^Monostilifera/)
			&& ($species !~ /^Scenedesmaceae/)
			&& ($species !~ /^Tubulinea/)
			&& ($species !~ /^Sordariomycetes/)
			&& ($species !~ /^Foraminifera/)
			&& ($species !~ /^Psocoptera/)
			&& ($species !~ /^Siphonostomatoida/)
			&& ($species !~ /^Tetraphyllidea/)
			&& ($species !~ /^Eurotiomycetes/)
			&& ($species !~ /^Chlorophyceae/)
			&& ($species !~ /^Microascaceae/)
			&& ($species !~ /^Haptoria/)
			&& ($species !~ /^Pleosporales/)
			&& ($species !~ /^Phaeothamniophyceae/)
			&& ($species !~ /^Enteropneusta/)
			&& ($species !~ /^Tubulanidae/)
			&& ($species !~ /^Diptera/)
			&& ($species !~ /^Tricladida/)
			&& ($species !~ /^Amoebozoa/)
			&& ($species !~ /^Bicosoecida/)
			&& ($species !~ /^Ectocarpales/)
			&& ($species !~ /^Pauropoda/)
			&& ($species !~ /^Pelagophyceae/)
			&& ($species !~ /^Chrysophyceae/)
			&& ($species !~ /^Eustigmatophyceae/)
			&& ($species !~ /^Thraustochytriidae/)
			&& ($species !~ /^Dactylogyridae/)
			&& ($species !~ /^Lasiosphaeriaceae/)
			&& ($species !~ /^Tardigrada/)
			&& ($species !~ /^Bdelloidea/)
			&& ($species !~ /^Crinoidea/)
			&& ($species !~ /^Rhizaria/)
			&& ($species !~ /^Spumella-like/)
			&& ($species !~ /^Siphonophorae/)
			&& ($species !~ /^Coccidia/)
			&& ($species !~ /^Chytridiales/)
			&& ($species !~ /^Roussoellaceae/)) { # A real genus name should be capitalized, else skip this part

		my @species = split('_',$species); # do this to accomodate for genera missing from SILVA tax file ex. Zea
		my $genusFromSpecies = shift(@species);
		my $genusFromLineage = pop(@lineage);

		if ($genusFromSpecies eq $genusFromLineage) {
			$observed{"genus"} = $genusFromLineage;
		}
		else {
			$observed{"genus"} = $genusFromSpecies;
			my $scalar = scalar(@lineage);
			$lineage[$scalar] = $genusFromLineage;
		}
	}
	
	foreach $taxon (@lineage) {
		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
		@taxon = split(/;/,$fulltaxon);
		$taxon = pop(@taxon); # use this to add to %observed

		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
			$rank = $map{$fulltaxon};
			$observed{$rank} = $taxon;
		}
		$prevtaxon = $fulltaxon;
	}
	$prevtaxon="";
}

##########

sub parse_header2_pro {
	%observed=();
	$line =~ s/^>//g;
	$line =~ s/_/ /; # replace FIRST underscore with space
	@line = split(' ',$line);
	$acc = shift(@line);
	$lineage = join(' ',@line);
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;

	# parse the original lineage
	@lineage = split(/;/,$lineage);
	$species = pop(@lineage);
	$observed{"species"} = $species;

	# grab the most abundant lineage from problem species	
	$lineage = $fixproblem{$species};
	$lineage =~ s/\s+/_/g; #replace spaces with underscores
	$lineage =~ s/\(//g; #remove parentheses
	$lineage =~ s/\)//g;
	$lineage =~ s/\.//g; #remove periods
	$lineage =~ s/\</_/g;
	$lineage =~ s/\>/_/g;
	@lineage = split(/;/, $lineage); # domain -> genus

	if (($species =~ /^[A-Z]/) && ($species !~ /^Eukaryota/) 
			&& ($species !~ /^Arthropoda/) 
			&& ($species !~ /^Tylenchida/)
			&& ($species !~ /^Cercozoa/)
			&& ($species !~ /^Polychaeta/)
			&& ($species !~ /^Ptychoderidae/)
			&& ($species !~ /^Palaeonemertea/)
			&& ($species !~ /^Trebouxiophyceae/)
			&& ($species !~ /^Chlorellales/)
			&& ($species !~ /^Chlorellaceae/)
			&& ($species !~ /^Peridiniales/)
			&& ($species !~ /^Harrimaniidae/)
			&& ($species !~ /^Ulvophyceae/)
			&& ($species !~ /^Dinobryon/)
			&& ($species !~ /^Heteronemertea/)
			&& ($species !~ /^Reticulitermes/)
			&& ($species !~ /^Cryptotermes/)
			&& ($species !~ /^Neotermes/)
			&& ($species !~ /^Myxozoa/)
			&& ($species !~ /^Solemya/)
			&& ($species !~ /^endosymbiont/)
			&& ($species !~ /^Rhinebothriidea/)
			&& ($species !~ /^Jakobida/)
			&& ($species !~ /^Nematoda/)
			&& ($species !~ /^Xanthophyceae/)
			&& ($species !~ /^Monostilifera/)
			&& ($species !~ /^Scenedesmaceae/)
			&& ($species !~ /^Tubulinea/)
			&& ($species !~ /^Sordariomycetes/)
			&& ($species !~ /^Foraminifera/)
			&& ($species !~ /^Psocoptera/)
			&& ($species !~ /^Siphonostomatoida/)
			&& ($species !~ /^Tetraphyllidea/)
			&& ($species !~ /^Eurotiomycetes/)
			&& ($species !~ /^Chlorophyceae/)
			&& ($species !~ /^Microascaceae/)
			&& ($species !~ /^Haptoria/)
			&& ($species !~ /^Pleosporales/)
			&& ($species !~ /^Phaeothamniophyceae/)
			&& ($species !~ /^Enteropneusta/)
			&& ($species !~ /^Tubulanidae/)
			&& ($species !~ /^Diptera/)
			&& ($species !~ /^Tricladida/)
			&& ($species !~ /^Amoebozoa/)
			&& ($species !~ /^Bicosoecida/)
			&& ($species !~ /^Ectocarpales/)
			&& ($species !~ /^Pauropoda/)
			&& ($species !~ /^Pelagophyceae/)
			&& ($species !~ /^Chrysophyceae/)
			&& ($species !~ /^Eustigmatophyceae/)
			&& ($species !~ /^Thraustochytriidae/)
			&& ($species !~ /^Dactylogyridae/)
			&& ($species !~ /^Lasiosphaeriaceae/)
			&& ($species !~ /^Tardigrada/)
			&& ($species !~ /^Bdelloidea/)
			&& ($species !~ /^Crinoidea/)
			&& ($species !~ /^Rhizaria/)
			&& ($species !~ /^Spumella-like/)
			&& ($species !~ /^Siphonophorae/)
			&& ($species !~ /^Coccidia/)
			&& ($species !~ /^Chytridiales/)
			&& ($species !~ /^Roussoellaceae/)) { # a real species would begin with a capital, else skip this part
		
		my @species = split('_',$species); # do this to accomodate for genera missing from SILVA tax file ex. Zea
		my $genusFromSpecies = shift(@species);
		my $genusFromLineage = pop(@lineage);

		if ($genusFromSpecies eq $genusFromLineage) {
			$observed{"genus"} = $genusFromLineage;
		}
		else {
			$observed{"genus"} = $genusFromSpecies;
			my $scalar = scalar(@lineage);
			$lineage[$scalar] = $genusFromLineage;
		}
	}

	# parse the most abundant lineage
	foreach $taxon (@lineage) {
		$fulltaxon = $prevtaxon.$taxon.";"; # use this to search %map
		@taxon = split(/;/,$fulltaxon);
		$taxon = pop(@taxon); # use this to add to %observed

		if (exists $map{$fulltaxon}) { #grab rank from map ( map only covers domain->genus )
			$rank = $map{$fulltaxon};
			$observed{$rank} = $taxon;
		}
		$prevtaxon = $fulltaxon;
	}

	$prevtaxon="";
}


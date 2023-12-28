#!/usr/bin/perl
#March 7, 2013 by Terri Porter
#script to replace use of ebot_nucleotide.plx and filter_results_numberingfixed.plx and filter_results_pickup.plx
#will need to add another script to do something similar to quick_filter_gb_results.plx
#Written to catch errors and move on instead of crashing the whole f-ing job
#USAGE perl test.plx

use strict;
use warnings;
use Bio::DB::EUtilities;

my $factory;

$factory = Bio::DB::EUtilities -> new (-eutil => 'esearch',
										-email => 'terriblue2002@yahoo.com',
										-db => 'nucleotide',
										-term => '("cox1" OR "coxI" OR "CO1" OR "COI") AND ("Acalyptus carpini"[Organism] OR "Acanthoscelides fraterculus"[Organism] OR "Acidota crenata"[Organism] OR "Acidota quadrata"[Organism] OR "Aclypea opaca"[Organism] OR "Actoribatella punctata"[Organism] OR "Acylophorus pronus"[Organism] OR "Adalia bipunctata"[Organism] OR "Aegialia cylindrica"[Organism] OR "Aegialia lacustris"[Organism] OR "Aegialia terminatis"[Organism] OR "Aelia americana"[Organism] OR "Agabus affinis"[Organism] OR "Agabus antennatus"[Organism] OR "Agabus anthracinus"[Organism] OR "Agabus approximatus"[Organism] OR "Agabus arcticus"[Organism] OR "Agabus bifarius"[Organism] OR "Agabus clavicornis"[Organism] OR "Agabus colymbus"[Organism] OR "Agabus confinis"[Organism] OR "Agabus discolor"[Organism] OR "Agabus elongatus"[Organism] OR "Agabus erichsoni"[Organism] OR "Agabus infuscatus"[Organism] OR "Agabus moestus"[Organism] OR "Agabus opacus"[Organism] OR "Agabus semipunctatus"[Organism] OR "Agabus seriatus"[Organism] OR "Agabus serricornis"[Organism] OR "Agabus strigulosus"[Organism] OR "Agabus thomsoni"[Organism] OR "Agabus tristis"[Organism] OR "Agathidium angulare"[Organism] OR "Agathidium jasperanum"[Organism] OR "Agonum affine"[Organism] OR "Agonum bicolor"[Organism] OR "Agonum consimile"[Organism] OR "Agonum corvus"[Organism] OR "Agonum cupreum"[Organism] OR "Agonum cupripenne"[Organism] OR "Agonum exaratum"[Organism] OR "Agonum gratiosum"[Organism] OR "Agonum harrisi"[Organism] OR "Agonum mutatum"[Organism] OR "Agonum propinquum"[Organism] OR "Agonum quinquepunctatum"[Organism] OR "Agonum retractum"[Organism] OR "Agonum simile"[Organism] OR "Agonum sordens"[Organism] OR "Agonum superioris"[Organism] OR "Agonum thoreyi"[Organism] OR "Agrotiphila staudingeri"[Organism] OR "Agrypnia straminea"[Organism] OR "Alaocybites egorovi"[Organism] OR "Aleochara castaneipennis"[Organism] OR "Aleochara suffusa"[Organism] OR "Aleochara verna"[Organism] OR "Altica ambiens"[Organism] OR "Altica bimarginata"[Organism] OR "Altica tombacina"[Organism] OR "Amara alpina"[Organism] OR "Amara bokori"[Organism] OR "Amara brunnea"[Organism] OR "Amara carinata"[Organism] OR "Amara colvillensis"[Organism] OR "Amara confusa"[Organism] OR "Amara discors"[Organism] OR "Amara ellipsis"[Organism] OR "Amara erratica"[Organism] OR "Amara farcta"[Organism] OR "Amara glacialis"[Organism] OR "Amara hyperborea"[Organism] OR "Amara interstitialis"[Organism] OR "Amara lacustris"[Organism] OR "Amara laevipennis"[Organism] OR "Amara littoralis"[Organism] OR "Amara obesa"[Organism] OR "Amara patruelis"[Organism] OR "Amara patruelis"[Organism] OR "Amara quelseli"[Organism] OR "Amara sinuosa"[Organism] OR "Amara spuria"[Organism] OR "Anachiptera latiteca"[Organism] OR "Anisotoma errans"[Organism] OR "Anthicus nigritus"[Organism] OR "Anthobium fimetarium"[Organism] OR "Aphodius albertanus"[Organism] OR "Aphodius borealis"[Organism] OR "Aphodius consentaneus"[Organism] OR "Aphodius corruptor"[Organism] OR "Aphodius guttatus"[Organism] OR "Aphodius pectoralis"[Organism] OR "Aphodius pectoralis"[Organism] OR "Aphodius piceatus"[Organism] OR "Aphodius tenellus"[Organism] OR "Aphodius yukonensis"[Organism] OR "Arctopora pulchella"[Organism] OR "Arctopsyche ladogensis"[Organism] OR "Arpedium cribratum"[Organism] OR "Asaphidion alaskanum"[Organism] OR "Asaphidion yukonense"[Organism] OR "Asaphidion yukonense"[Organism] OR "Asaphidion yukonense"[Organism] OR "Atheta granulata"[Organism] OR "Atheta hyperborea"[Organism] OR "Atomaria kamtschatica"[Organism] OR "Atomaria pusilla"[Organism] OR "Auleutes epilobii"[Organism] OR "Auleutes nebulosus"[Organism] OR "Bembidion aquiliferum"[Organism] OR "Bembidion arcticum"[Organism] OR "Bembidion balli"[Organism] OR "Bembidion bimaculatum"[Organism] OR "Bembidion breve"[Organism] OR "Bembidion chalceum"[Organism] OR "Bembidion coloradense"[Organism] OR "Bembidion complanulum"[Organism] OR "Bembidion compressum"[Organism] OR "Bembidion concolor"[Organism] OR "Bembidion concretum"[Organism] OR "Bembidion convexulum"[Organism] OR "Bembidion curtulatum"[Organism] OR "Bembidion dauricum"[Organism] OR "Bembidion dyschirinum"[Organism] OR "Bembidion dyschirinum"[Organism] OR "Bembidion fortestriatum"[Organism] OR "Bembidion foveum"[Organism] OR "Bembidion grapii"[Organism] OR "Bembidion grapii"[Organism] OR "Bembidion hasti"[Organism] OR "Bembidion honestum"[Organism] OR "Bembidion hyperboraeorum"[Organism] OR "Bembidion insulatum"[Organism] OR "Bembidion intermedium"[Organism] OR "Bembidion interventor"[Organism] OR "Bembidion lapponicum"[Organism] OR "Bembidion levettei"[Organism] OR "Bembidion morulum"[Organism] OR "Bembidion mutatum"[Organism] OR "Bembidion nigripes"[Organism] OR "Bembidion nigrum"[Organism] OR "Bembidion nitidum"[Organism] OR "Bembidion patruele"[Organism] OR "Bembidion petrosum"[Organism] OR "Bembidion planatum"[Organism] OR "Bembidion poppi"[Organism] OR "Bembidion pseudocautum"[Organism] OR "Bembidion quadrifoveolatum"[Organism] OR "Bembidion quadrimaculatum"[Organism] OR "Bembidion roosvelti"[Organism] OR "Bembidion salebratum"[Organism] OR "Bembidion scudderi"[Organism] OR "Bembidion sejunctum"[Organism] OR "Bembidion semipunctatum"[Organism] OR "Bembidion semistriatum"[Organism] OR "Bembidion sordidum"[Organism] OR "Bembidion sordidum"[Organism] OR "Bembidion transparens"[Organism] OR "Bembidion umiatense"[Organism] OR "Bembidion yukonum"[Organism] OR "Berninelsonius hyperboreus"[Organism] OR "Bledius confusus"[Organism] OR "Bledius viriosus"[Organism] OR "Bledius zophus"[Organism] OR "Blethisa catenaria"[Organism] OR "Blethisa inexspectata"[Organism] OR "Blethisa multipunctata"[Organism] OR "Bolitobius analis"[Organism] OR "Boreaphilus henningianus"[Organism] OR "Boreophilia fusca"[Organism] OR "Boreophilia gelida"[Organism] OR "Boreophilia nearctica"[Organism] OR "Bradycellus lecontei"[Organism] OR "Bromius obscurus"[Organism] OR "Byrrhus cyclophorus"[Organism] OR "Byrrhus eximius"[Organism] OR "Byrrhus fasciatus"[Organism] OR "Byrrhus kirbyi"[Organism] OR "Byturus unicolor"[Organism] OR "Caenocara scymnoides"[Organism] OR "Calligrapha californica"[Organism] OR "Camponotus herculeanus"[Organism] OR "Cantharis aneba"[Organism] OR "Carabus chamissonis"[Organism] OR "Carabus maeander"[Organism] OR "Carabus taedatus"[Organism] OR "Carabus truncaticollis"[Organism] OR "Carabus truncaticollis"[Organism] OR "Carabus vietinghoffi"[Organism] OR "Carphoborus andersoni"[Organism] OR "Carphoborus carri"[Organism] OR "Cassida flaveola"[Organism] OR "Catops alpinus"[Organism] OR "Catops alsiosus"[Organism] OR "Catops americanus"[Organism] OR "Catops basilaris"[Organism] OR "Cepheus corae"[Organism] OR "Ceratomegilla ulkei"[Organism] OR "Ceratoppia bipilis"[Organism])',
										-usehistory => 'y');

my $count = $factory -> get_count;

#get history from queue
my $hist = $factory -> next_History || die 'No history data returned';
print "History returned\n";

#db carries over from above
$factory -> set_parameters (-eutil => 'efetch',
							-rettype => 'gb',
							-history => $hist);

my $retry = 0;
my ($retmax, $retstart) = (500,0);

open (my $out, ">", "seqs.gb") || die "Can't open file: $!\n";

RETRIEVE_SEQS:
while ($retstart < $count) {
	$factory -> set_parameters (-retmax => $retmax,
								-retstart => $retstart);
	
	eval {
		$factory -> get_Response (-cb => sub {
											my ($data) = @_;
											print $out $data
										 }
								 );
	};

	if ($@) {
		die "Server error : $@.  Try again later" if $retry == 5;
		print STDERR "Server error, redo #$retry\n";
		$retry++ && redo RETRIEVE_SEQS;
	}

	print "Retrieved $retstart\n";
	$retstart += $retmax;

}

close $out;


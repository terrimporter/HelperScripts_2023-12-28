#!/usr/bin/perl
#Feb. 22, 2012 numbering of gb outfiles fixed to correspond with batch number
#March 27, 2012 by Terri Porter
#This script will grab reference sequences, more than one per taxon!
#used ebot to grab list of taxids and gis
#ebot_taxonomy.plx "Vertebrata"[Organism] AND "species"[rank] > taxonomy.taxid
#ebot_nucleotide.plx "cytochrome oxidase subunit I" OR "cytochrome c oxidase subunit I" > nucleotide.gi
#usage perl filter_results.plx taxonomy.taxid nucleotide.gi

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

#declare var
my $num_gi;
my $num_parts;
my $num_parts_roundup;
my $i=1;
my $name;
my $item;
my $j=0;
my $file;
my $gi_batch;
my $string="";

#declare array
my @gi;
my @gi_part;
my @batch_files;
my @gi_batch;

open (GI,"<",$ARGV[1]) || die "Error reading nucleotide.gi: $!\n";
@gi = <GI>;
close GI;

$num_gi = scalar(@gi);
$num_parts = $num_gi/500;#edit batch size here
$num_parts_roundup = $num_parts + 1;

#split gi list into batches of 500 gi numbers each
while ($i <= $num_parts_roundup) {
	@gi_part = splice(@gi,0,500); #edit batch size here too
	$name = "gi_".$i.".batch";
	open (OUT,">>",$name) || die "Cannot write to $name: $!\n";

	while ($gi_part[$j]) {
		$item = $gi_part[$j];
		chomp $item;
		print OUT "$item\n";
		$j++;
	}
	$i++;
	$j=0;
	@gi_part=();
	close OUT;
}
$i=0;

#put all batch files into array
@batch_files = qx(ls | grep '.batch');

while ($batch_files[$i]) {
	$file = $batch_files[$i];
	chomp $file;

	open (BATCH,"<",$file) || die "Canot open $file: $!\n";
	@gi_batch = <BATCH>;
#	print "gi_batch array @gi_batch\n";#test

	while ($gi_batch[$j]) {
		$gi_batch = $gi_batch[$j];
		chomp $gi_batch;
		if ($j==0) {
			$string = $gi_batch;
		}
		else {
			$string = $string.",".$gi_batch;
		}
		$j++;
	}
	$j=0;
	$gi_batch=();

#	print "$string";#test

	run_efetch();
	
	close BATCH;
	$i++;
	@gi_batch=();
	$string="";
}
$i=0;


####################

sub run_efetch {

#declare var
my $taxid;
my $factory;
#my $file;
my $num;
my $file2;
my $seqin;
my $seq;
my $gb;
my $feat_object;
my $value;
my $check;
my $organism;
my $sequence_fragment;
my $COI_seq;
my $gi;

#declare array
my @taxid;
my @ids;

#declare hash
my %taxid;

open (TAXID,"<",$ARGV[0]) || die "Error reading taxonomy.taxid: $!\n";
@taxid = <TAXID>;
close TAXID;

while ($taxid[$j]) {
	$taxid = $taxid[$j];
	chomp $taxid;
	$taxid{$taxid} = 1;
	$j++;
}
$j=0;

#while( my($key,$value) = each(%taxid)) { #test
#	print "$key	=>	$value\n";
#}

$factory = Bio::DB::EUtilities -> new (	-eutil	=>	'efetch',
										-db		=>	'nucleotide',
										-rettype	=>	'gb',
										-retmode	=> 'text',
										-email	=> 'terriblue2002@yahoo.com',
										-id	=>	$string);#submit as string instead of array

#$file2 = 'cox1_'.$i.'.gb';
$file =~ /gi_(\d+)\.batch/;
$num = $1;
$file2 = 'cox1_'.$num.'.gb';

$factory -> get_Response (-file	=>	$file2);

$seqin = Bio::SeqIO -> new(	-file	=> $file2,
							-format	=> 'genbank');

open (GB_TAXID,">>","gb_taxid.map") || die "Error cannot write to gb_taxid.map: $!\n";
print GB_TAXID "GB\tTaxId\n";

open (GB_ORG,">>","gb_org.map") || die "Error cannot write to gb_org.map: $!\n";
print GB_ORG "GB\tOrganism\n";


open (GB_SEQ,">>","gb_seq.map") || die "Error cannot write to gb_seq.map: $!\n";
print GB_SEQ "GB\tSequence\n";

open (GB_GI,">>","gb_gi.map") || die "Error cannot write to gb_gi.map: $!\n";
print GB_GI "GB\tGI\n";

while ($seq = $seqin -> next_seq) {
	$gb = $seq -> id;
	for $feat_object ($seq -> get_SeqFeatures) {
		if ($feat_object -> primary_tag eq "source") {
			if ($feat_object -> has_tag('db_xref')) {
				for $value ($feat_object -> get_tag_values('db_xref')) {
					if ($value =~ /taxon:/) {
						$value =~ s/taxon://;
						$taxid = $value;
						print GB_TAXID "$gb\t$taxid\n";
						#print "taxid $value\t";#test
						$check = $taxid{$taxid}; #filter by taxonomy
						#print "check $check\n";#test
						if ($check == 1) {
							#$taxid{$taxid}=0;#only retrieve one sequence per taxon
							if ($feat_object -> has_tag('organism')) {
								for $value ($feat_object -> get_tag_values('organism')) {
									$organism = $value;
									print GB_ORG "$gb\t$organism\n";
								}
							}
						}
					}
				}
			}
		}
		if ($check ==1) {  #change from elsif to if statement!!!! 12/23/12
			if ($feat_object -> primary_tag eq "CDS") {
				$sequence_fragment = $feat_object -> spliced_seq -> seq; #grab COI nt seq
				#print "sequence fragment $sequence_fragment\n";#test
				if ($feat_object -> has_tag('gene')) {
					for $value ($feat_object -> get_tag_values('gene')) {
						if ($value =~ /(cox1|coxI|CO1|COI)/i) { #double check product ###added CO1 and COI
						#print "value gene $value\n";#test
						$COI_seq = $sequence_fragment;
						#$taxid{$taxid}=0;#only retrieve one sequence per taxon
						print GB_SEQ "$gb\t$COI_seq\n";
						#print FEATURES "$id\t$taxid\t$organism\t$COI_seq\n";
						}
					}
				}
				if ($feat_object -> has_tag('db_xref')) {
					for $value ($feat_object -> get_tag_values('db_xref')) {
						if ($value =~ /GI:/) {
							$value =~ s/GI://;
							$gi = $value;
							print GB_GI "$gb\t$gi\n";
						}
					}
				}

			}
			#print FEATURES "$id\t$taxid\t$organism\t$COI_seq\n";
		}
	}
	#print FEATURES "$id\t$gi\t$taxid\t$organism\t$COI_seq\n";
	$gb=();
	$feat_object=();
	$value=();
	$taxid=();
	$check=();
	$organism=();
	$sequence_fragment=();
	$COI_seq=();
}
close GB_TAXID;
close GB_ORG;
close GB_SEQ;
close GB_GI;

}

####################

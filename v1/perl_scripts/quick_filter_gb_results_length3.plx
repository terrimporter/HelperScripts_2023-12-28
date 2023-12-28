#!/usr/bin/perl
#May 9, 2012 updated to filter through a list of existing gb files in dir
#NEW USAGE perl quick_filter_results.plx taxonomy.taxid

#May 2, 2012 updated .gi and .batch numbers to be equivalent
#March 27, 2012 by Terri Porter
#used ebot to grab list of taxids and gis
#ebot_taxonomy.plx "Vertebrata"[Organism] AND "species"[rank] > taxonomy.taxid
#ebot_nucleotide.plx "cytochrome oxidase subunit I" OR "cytochrome c oxidase subunit I" > nucleotide.gi
#usage perl filter_results.plx taxonomy.taxid nucleotide.gi

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

#declare var
my $taxid;
my $j=0;
my $dir;
my $i=0;
my $file;
my $filepath;
my $num;
my $file2;
my $factory;
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
my $length;
my $minlength=500; #### reset this if necessary to 250 for trnL, otherwise 500bp) #####

#declare array
my @taxid;
my @files;
my @COI_seq;

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

print "Enter path to dir containing gb files including final /:\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error opening dir $dir\n";
@files = readdir DIR;
closedir DIR;

while ($files[$i]) {
	$file = $files[$i];

	if ($file =~ /^\./) {
		$i++;
		next;
	}
	else {

		$filepath = $dir.$file;

		$seqin = Bio::SeqIO -> new(	-file	=> $filepath,
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
								if (exists $taxid{$taxid}) {
									#$check = $taxid{$taxid}; #filter by taxonomy
									#print "check $check\n";#test
									#if ($check == 1) {
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
				if ($feat_object -> primary_tag eq "CDS") {
					$sequence_fragment = $feat_object -> spliced_seq -> seq; #grab COI nt seq
					#print "sequence fragment $sequence_fragment\n";#test
					if ($feat_object -> has_tag('gene')) {
						for $value ($feat_object -> get_tag_values('gene')) {
							if ($value =~ /(cox1|coxI|CO1|COI)/i) { #double check product
								#print "value gene $value\n";#test
								$COI_seq = $sequence_fragment;
								@COI_seq = split(//,$COI_seq);
								$length = scalar(@COI_seq);
								if ($length >= $minlength) {
									#$taxid{$taxid}=0;#only retrieve one sequence per taxon
									print GB_SEQ "$gb\t$COI_seq\n";
								}
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
			}
			$gb=();
			$feat_object=();
			$value=();
			$taxid=();
			$check=();
			$organism=();
			$sequence_fragment=();
			$COI_seq=();
			$seq=();
			$gi=();
		}
		$num=();
		$file2=();
		$factory=();
		$seqin=();

	}
	$i++;
	$file=();
}
$i=0;
close GB_TAXID;
close GB_ORG;
close GB_SEQ;
close GB_GI;

####################

#!/usr/bin/perl
#April 10, 2018 clean up this script a bit, edit to produce outfiles with modified filenames so they don't get mixed up with the original mapping files, also adjust gene regex to find bacterial ctaD 
#USAGE perl quick_filter_gb_results_length5.plx

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
my $value_dbxref;
my $value_organism;
my $value_gene;
my $value_dbxref2;

#declare array
my @taxid;
my @files;
my @COI_seq;
my @value_dbxref;
my @feat_object;
my @value_organism;
my @value_gene;
my @value_dbxref2;

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

print "Enter path to dir containing gb files including final /:\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error opening dir $dir\n";
@files = readdir DIR;
closedir DIR;

open (GB_TAXID,">>","gb_taxid_bact.map") || die "Error cannot write to gb_taxid_bact.map: $!\n";

open (GB_SEQ,">>","gb_seq_bact.map") || die "Error cannot write to gb_seq_bact.map: $!\n";


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

#		open (GB_TAXID,">>","gb_taxid_bact.map") || die "Error cannot write to gb_taxid_bact.map: $!\n";
#		print GB_TAXID "GB\tTaxId\n";

#		open (GB_SEQ,">>","gb_seq_bact.map") || die "Error cannot write to gb_seq_bact.map: $!\n";
#		print GB_SEQ "GB\tSequence\n";

		while ($seq = $seqin -> next_seq) {
			$gb = $seq -> id;
			@feat_object = $seq -> get_SeqFeatures;
			foreach $feat_object (@feat_object) {
				if ($feat_object -> primary_tag eq "source") {
					if ($feat_object -> has_tag('db_xref')) {
						@value_dbxref = $feat_object -> get_tag_values('db_xref'); 
						foreach $value_dbxref (@value_dbxref) {
							if ($value_dbxref =~ /taxon:/) {
								$value_dbxref =~ s/taxon://;
								$taxid = $value_dbxref;
								print GB_TAXID "$gb\t$taxid\n";
							}
						}
					}
				}
				if ($feat_object -> primary_tag eq "CDS") {
					$sequence_fragment = $feat_object -> spliced_seq -> seq; #grab COI nt seq
					if ($feat_object -> has_tag('gene')) {
						@value_gene = $feat_object -> get_tag_values('gene');
						foreach $value_gene (@value_gene) {
							if ($value_gene =~ /^ctaD$/i) { #double check product
								$COI_seq = $sequence_fragment;
								@COI_seq = split(//,$COI_seq);
								$length = scalar(@COI_seq);
								if ($length >= $minlength) {
									print GB_SEQ "$gb\t$COI_seq\n";
								}
							}
						}
					}
				}
			}
			$gb=();
			$feat_object=();
			$value_dbxref=();
			$value_gene=();
			$taxid=();
			$check=();
			$sequence_fragment=();
			$COI_seq=();
			$seq=();
			@value_dbxref=();
			@feat_object=();
			@value_gene=();
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
close GB_SEQ;

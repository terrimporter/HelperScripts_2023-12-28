#!/usr/bin/perl
#Sept. 16, 2011 edit to work with ITS gi's from ITS_entrez.plx
#Sept. 15, 2011 use gi.list and grab genbank records in batches of 500
#usage perl ITS_entrez_grab_gb.plx gi.query

use strict;
use warnings;
use Bio::SeqIO;
use Bio::DB::EUtilities;

#declare var
my $i=1;
my $num_gi;
my $num_parts;
my $num_parts_roundup;
my $j=0;
my $k=0;
my $file;
my $l=0;
my $name;
my $item;

#declare array
my @gi;
my @gi_part;
my @batch_files;
my @gi_batch;

#read in gi list from ITS_entrez.plx

open (IN,"<",$ARGV[0]) || die ("Error cannot read infile: $!\n");
@gi = <IN>;
close IN;

$num_gi = scalar(@gi);
$num_parts = $num_gi/500; #edit batch size here
$num_parts_roundup = $num_parts + 1;

#split list into batches of 500 gi numbers each
while ($i <= $num_parts_roundup) {
	@gi_part = splice(@gi,0,500); #edit batch size here too
	$name = "gi_".$i.".batch";
	open (OUT,">>",$name) || die "Can't open $name: $!\n";

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

#put all batch files into an array to process one at a time
@batch_files = qx(ls | grep '.batch');

while ($batch_files[$k]) {
	$file = $batch_files[$k];
	chomp $file;

	open (BATCH,"<",$file) || die "Can't open $file: $!\n";
	@gi_batch = <BATCH>;

	run_efetch();
	
	close BATCH;
	$k++;
}

####################

sub run_efetch {

#declare var
my $factory;
my $file2;
my $seqin;
my $seq;
my $id;
my $feat_object;
my $organism;
my $sequence_fragment;
my $value;
my $ITS_seq="nil";
my $spacer1="nil";
my $spacer2="nil";

#declare array

$factory = Bio::DB::EUtilities -> new (	-eutil	=>	'efetch',
					-db	=>	'nucleotide',
					-rettype=>	'gb',
					-id	=>	\@gi_batch);
$file2 = 'ITS_'.$k.'.genbank';

$factory -> get_Response (	-file	=> $file2);
$seqin = Bio::SeqIO -> new (	-file	=> $file2,
				-format	=> 'genbank');

open (FEATURES,">>","features.txt") || die ("Error cannot write to features.txt: $!\n");
print FEATURES "ID\tOrganism\tITS sequence\n";

while ($seq = $seqin -> next_seq) {
	$id = $seq -> id;
	for $feat_object ($seq -> get_SeqFeatures) {
		if ($feat_object -> primary_tag eq "source") {
			if ($feat_object -> has_tag('organism')) {
				for $value ($feat_object ->get_tag_values('organism')) {
					$organism = $value;
				}
			}
		}
		elsif ($feat_object -> primary_tag eq "rRNA") {
			$sequence_fragment = $feat_object -> spliced_seq -> seq;
			if ($feat_object -> has_tag('product')) {
				for $value ($feat_object -> get_tag_values('product')) {
					if ($value =~ /5.8S/) {
						$ITS_seq = $sequence_fragment;
					}
				}
			}
		}
		elsif ($feat_object->primary_tag eq "misc_RNA") {
			$sequence_fragment = $feat_object -> spliced_seq->seq;
			if ($feat_object->has_tag('product')) {
				for $value ($feat_object->get_tag_values('product')) {
					if ($value =~/spacer 1/) {
						$spacer1 = $sequence_fragment;
					}
					elsif ($value =~ /spacer 2/) {
						$spacer2 = $sequence_fragment;
					}
				}
			}
		}
	}
	print FEATURES "$id\t$organism\t$spacer1\t$ITS_seq\t$spacer2\n";
	$ITS_seq="nil";
	$spacer1="nil";
	$spacer2="nil";
}
close FEATURES;

}

####################

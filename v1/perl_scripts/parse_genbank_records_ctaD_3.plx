#!/usr/bin/perl
#July 4, 2020 based on quick_filter_gb_results_length5b.plx
# Script to parse through genbank records using Bio::SeqIO
#USAGE perl parse_genbank_records_ctaD.plx

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Number::Range;

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
my $range;
my $lower_numeric;
my $upper_numeric;
my $outfile;

#declare array
my @taxid;
my @files;
my @COI_seq;
my @value_dbxref;
my @feat_object;
my @value_organism;
my @value_gene;
my @value_dbxref2;
my @range;

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

		open (GB_TAXID,">>","gb_taxid_ctaD.map") || die "Error cannot write to gb_taxid.map: $!\n";

		open (GB_SEQ,">>","gb_seq_ctaD.map") || die "Error cannot write to gb_seq.map: $!\n";

		while ($seq = $seqin -> next_seq) {
			$gb = $seq -> id;
			# parse feature table when available 
			@feat_object = $seq -> get_SeqFeatures;

			parse_seq_features(\@feat_object);

			# if no seq was processed above, then parse annotations
			if (! defined $COI_seq) {
				my $anno_collection = $seq -> annotation;
				for my $key ( $anno_collection -> get_all_annotation_keys ) {
    				my @annotations = $anno_collection->get_Annotations('wgs');
    				foreach my $value ( @annotations ) {
						$range= $value -> display_text;
						@range=split(/-/,$range); #reformat for perl range
						my $lower = $range[0];
						my $upper = $range[1];
						if ($lower =~ /^([A-Z]+)([0-9]+)$/ ) {
							my $lower_prefix = $1;
							my $lower_numeric = $2;
							
							if ($upper =~ /^([A-Z]+)([0-9]+)$/) {
								my $upper_prefix = $1;
								my $upper_numeric = $2;
								my @arr2;
								for ($lower_numeric..$upper_numeric) {
									my $newval = $upper_prefix.$_;
									push(@arr2,$newval);
								}

								# grab all wgs records
								$outfile = submit_entrez_fetch(\@arr2,$range);

								# parse features from all wgs records
								my $seqin2 = Bio::SeqIO -> new(	-file	=> $outfile,
															-format	=> 'genbank');

								while ($seq = $seqin2 -> next_seq) {
									$gb = $seq -> id;
									# parse feature table when available 
									@feat_object = $seq -> get_SeqFeatures;

									parse_seq_features(\@feat_object);

								}
							}
						}
    				}
				}
				$range=();
				@range=();
				unlink $outfile; # these files are huge, get rid of them after parsing
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

#####

sub parse_seq_features {

	@feat_object = @{$_[0]};

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
			if ($feat_object -> has_tag('organism')) {
				$organism = $feat_object -> get_tag_values('organism');
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
						if ($length >= $minlength) { #double check length
							print GB_SEQ "$gb\t$COI_seq\n";
						}
					}
				}
			}
		}
	}
}

#####

sub submit_entrez_fetch {


	my @arr2 = @{$_[0]};
	my $range = $_[1];

	$outfile = $range.".gb";

	my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
      	                                  -db      => 'nucleotide',
           		                           -id      => \@arr2,
                   		                   -email   => 'terriblue2002@yahoo.com',
                           		           -rettype => 'gbwithparts');

	# file with sequence
	$factory->get_Response(-file => $outfile);

	return($outfile);

}






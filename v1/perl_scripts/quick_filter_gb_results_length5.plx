#!/usr/bin/perl
#Sept. 16, 2016; Rewrite script to parse one GenBank file at a time, use parallel to pass filenames, and output errors
#NEWUSAGE: ls | grep seqs.gb | parallel -j 1 "perl quick_filter_gb_results_length5g.plx taxonomy.taxid {}"
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
my $file;
my $seqin=();
my $seq;
my $gb;
my $feat_object;
my $value_organism;
my $organism;
my $value_dbxref;
my $sequence_fragment;
my $value_gene;
my $COI_seq;
my $length;
my $minlength=500; #### reset this if necessary to 250 for trnL, otherwise 500bp) #####
my $value_dbxref2;
my $gi;

#declare array
my @taxid;
my @feat_object;
my @value_dbxref;
my @value_organism;
my @value_gene;
my @COI_seq;
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

open (LOG, ">>", "log.txt") || die "Error cannot write to log.txt: $!\n";

open (GB_TAXID,">>","gb_taxid.mapredo") || die "Error cannot write to gb_taxid.map: $!\n";
print GB_TAXID "GB\tTaxId\n";

open (GB_ORG,">>","gb_org.mapredo") || die "Error cannot write to gb_org.map: $!\n";
print GB_ORG "GB\tOrganism\n";

open (GB_SEQ,">>","gb_seq.mapredo") || die "Error cannot write to gb_seq.map: $!\n";
print GB_SEQ "GB\tSequence\n";

open (GB_GI,">>","gb_gi.mapredo") || die "Error cannot write to gb_gi.map: $!\n";
print GB_GI "GB\tGI\n";

$file = $ARGV[1];
chomp $file;

if ($file eq '.' || $file eq '..') { #use a different method to skip . and ..
print LOG "Skipping file: ".$file."\n";#test
}
else {
	$seqin = Bio::SeqIO -> new(	-file	=> $file,
								-format	=> 'genbank');

	if (defined $seqin) {#no warning if $seqin is defined
		print LOG "\tProcessing file: ".$file."\n"; #test
	
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
									if (exists $taxid{$taxid}) {
										if ($feat_object -> has_tag('organism')) {
											@value_organism = $feat_object -> get_tag_values('organism');
											foreach $value_organism (@value_organism) {
												$organism = $value_organism;
												print GB_ORG "$gb\t$organism\n";
												$value_organism=();
												$organism=();
											}
											@value_organism=();
										}
									}
									$taxid=();
								}
								$value_dbxref=();
							}
							@value_dbxref=();
						}
					}
					if ($feat_object -> primary_tag eq "CDS") {				
					$sequence_fragment = $feat_object -> spliced_seq -> seq; #grab COI nt seq
					#print "sequence fragment $sequence_fragment\n";#test
						if ($feat_object -> has_tag('gene')) {
							@value_gene = $feat_object -> get_tag_values('gene');
							foreach $value_gene (@value_gene) {
								if ($value_gene =~ /(cox1|coxI|CO1|COI)/i) { #double check product
									#print "value gene $value\n";#test
									#$sequence_fragment = $feat_object -> spliced_seq -> seq; #grab COI nt seq
									$COI_seq = $sequence_fragment;
									@COI_seq = split(//,$COI_seq);
									$length = scalar(@COI_seq);
									if ($length >= $minlength) {
										#$taxid{$taxid}=0;#only retrieve one sequence per taxon
										print GB_SEQ "$gb\t$COI_seq\n";
									}
									#print FEATURES "$id\t$taxid\t$organism\t$COI_seq\n";
									$COI_seq=();
									@COI_seq=();
									$length=();
									#$sequence_fragment=();
								}
								$value_gene=();
							}
							@value_gene=();
						}
						if ($feat_object -> has_tag('db_xref')) {
							@value_dbxref2 = $feat_object -> get_tag_values('db_xref');
							foreach $value_dbxref2 (@value_dbxref2) {
								if ($value_dbxref2 =~ /GI:/) {
									$value_dbxref2 =~ s/GI://;
									$gi = $value_dbxref2;
									print GB_GI "$gb\t$gi\n";
								}
								$value_dbxref2=();
								$gi=();
							}
							@value_dbxref2=();
						}
						$sequence_fragment=();
					}
					$feat_object=();
				}
				$gb=();
				$seq=();
				@feat_object=();
			}
			
	}
	else {
		print "\t\tProblem processing file $file:\n"; #test
	}
	$seqin=();
}
close GB_TAXID;
close GB_ORG;
close GB_SEQ;
close GB_GI;
close LOG;

####################

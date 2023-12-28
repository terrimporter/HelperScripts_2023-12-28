#!/usr/bin/perl
#Dec. 8/17 get /country and /lat_lon for a LIST of directories in a .txt file including final '/' for each dir
#USAGE perl quick_filter_country_latlon_CO1.plx taxonomy.taxid dirlist.txt

#Nov. 7/17 updated to get FEATURE source /country for all CO1 sequences

#May 9, 2012 updated to filter through a list of existing gb files in dir
#NEW USAGE perl quick_filter_results.plx taxonomy.taxid

#May 2, 2012 updated .gi and .batch numbers to be equivalent
#March 27, 2012 by Terri Porter
#used ebot to grab list of taxids and gis
#ebot_taxonomy.plx "Vertebrata"[Organism] AND "species"[rank] > taxonomy.taxid
#ebot_nucleotide.plx "cytochrome oxidase subunit I" OR "cytochrome c oxidase subunit I" > nucleotide.gi
#usage perl quick_filter_gb_country_latlon_CO1.plx taxonomy.taxid

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

#declare var
my $taxid;
my $j=0;
my $dir;
my $year;
my $i=0;
my $file;
my $filepath;
my $outfilename;
my $seqin;
my $seq;
my $gb;
my $feat_object;
my $country="nil";
my $latlon="nil";
my $sequence_fragment;
my $COI_seq;
my $gi;
my $length;
my $minlength=500; #### reset this if necessary to 250 for trnL, otherwise 500bp) #####
my $value_gene;
my $value;
my $name;

#declare array
my @dirlist;
my @dirarray;
my @taxid;
my @files;
my @COI_seq;
my @value_country;
my @feat_object;
my @country;
my @value_gene;

#declare hash
my %taxid;
#my %hash; #key1=gb, key2=country, key3=CO1seq500bp+ (binary 0/1)

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

open (DIRLIST,"<",$ARGV[1]) || die "Error reading dirlist.txt: $!\n";
@dirlist = <DIRLIST>;
close DIRLIST;

foreach $dir (@dirlist) {
	chomp $dir;

	@dirarray = split('_',$dir);
	$year = $dirarray[2];
	$year =~ s/\///g;

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

			$outfilename = 'country_latlon_'.$year.'.txt';

			open (OUT,">>",$outfilename) || die "Error cannot write to $outfilename: $!\n";

			while ($seq = $seqin -> next_seq) {
				$gb = $seq -> id;
				@feat_object = $seq -> get_SeqFeatures;
				foreach $feat_object (@feat_object) {
					if ($feat_object -> primary_tag eq "source") {
						if ($feat_object -> has_tag('country')) {
							for $value ($feat_object -> get_tag_values('country')){
								$country = $value;
							}
						}
						else {
							$country="nil";
						}
						if ($feat_object -> has_tag('lat_lon')) {
							for $value ($feat_object -> get_tag_values('lat_lon')){
								$latlon = $value;
							}
						}
						else {
							$latlon="nil";
						}
					}
					if ($feat_object -> primary_tag eq "CDS") {
						$sequence_fragment = $feat_object -> spliced_seq -> seq; #grab COI nt seq
						if ($feat_object -> has_tag('gene')) {
							@value_gene = $feat_object -> get_tag_values('gene');
							foreach $value_gene (@value_gene) {
								if ($value_gene =~ /^(cox1|coxI|CO1|COI)$/i) { #double check product
									$COI_seq = $sequence_fragment;
									@COI_seq = split(//,$COI_seq);
									$length = scalar(@COI_seq);
									if ($length >= $minlength) { #only count if 500bp+
										if ($country eq "nil") {
											if ($latlon eq "nil") {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t1\t0\t0\n";
											}
											else {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t1\t0\t1\n";
											}
										}
										else {
											if ($latlon eq "nil") {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t1\t1\t0\n";
											}
											else {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t1\t1\t1\n";
											}
										}
									}
									else {
										if ($country eq "nil") {
											if ($latlon eq "nil") {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t0\t0\t0\n";
											}
											else {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t0\t0\t1\n";
											}
										}
										else {
											if ($latlon eq "nil") {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t0\t1\t0\n";
											}
											else {
												print OUT "species\t$gb\t$year\t$country\t$latlon\t0\t1\t1\n";
											}
										}
									}
								}
							}
						}
					}
				}
				$gb=();
				$country="nil";
				$latlon="nil";
				$feat_object=();
				$value_gene=();
				$taxid=();
				$sequence_fragment=();
				$COI_seq=();
				$seq=();
				@feat_object=();
				@value_gene=();
			}
			$seqin=();
			close OUT;
		}
		$i++;
		$file=();
	}
	$i=0;
	@files=();
}

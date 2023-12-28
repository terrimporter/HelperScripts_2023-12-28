#!/usr/bin/perl
# Teresita M. Porter, March 16, 2020
# Script to parse annotation and FASTA files from MitoFish 
# Original data from http://mitofish.aori.u-tokyo.ac.jp/download.html
# Annotation .txt files in mitogenome_annotations dir
# Mitogenome .fa files in all_mitogenomes dir
# USAGE perl parse_mitofish.plx mitogenome_annotations all_mitogenomes

use strict;
use warnings;
use Data::Dumper;

#var
my $annot_dir = $ARGV[0];
my $fa_dir = $ARGV[1];
my $i=0;
my $line;
my $id_line;
my $molec_type;
my $range;
my $gene_product;
my $id;
my $start;
my $stop;
my $outfile = "12S.fasta";
my $filename;
my $filepath;
my $j=0;
my $seq="";
my $length;
my $X12S;

#array
my @annot;
my @line;
my @id_line;
my @range;
my @in;
my @fa_filenames;
my @filename;


#hash
my %range; # key = id, value = range

# grab annotation info
chdir $annot_dir;
@annot = `grep "12S" *`;
chdir "..";

# parse annotations
while ($annot[$i]) {
	$line = $annot[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$id_line = $line[0];
	$molec_type = $line[1];
	$range = $line[2];
	$gene_product = $line[3];

	@id_line = split(/\./, $id_line);
	$id = $id_line[0];

	$range{$id} = $range;

	$i++;

}
$i=0;

# grab fa filenames
opendir (DIR2, $fa_dir) || die $!;
@fa_filenames = readdir (DIR2);
closedir (DIR2) || die $!;

open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

# open each FASTA file, get 12S region, print new 12S FASTA
while ($fa_filenames[$i]) {
	$filename = $fa_filenames[$i];
	chomp $filename;

	# skip the FASTA files that contain just the genes since not all taxa have this
	unless ($filename =~ /genes/ || $filename =~/^\./) {

		# get id from filename
		@filename = split(/\./,$filename);
		$id = $filename[0];

		$filepath = $fa_dir."/".$filename;
		open (IN, "<", $filepath) || die "Cannot open fasta file: $!\n";
		@in = <IN>;
		close IN;

		while ($in[$j]) {
			$line = $in[$j];
			chomp $line;

			unless ($line =~ /^>/) {
				$seq = $seq.$line;
			}
			$j++;
		}
		$j=0;

		if (exists $range{$id}) {
			$range = $range{$id};
			@range = split(/\.\./, $range);
			$start = $range[0];
			$start = $start-1;
			$stop = $range[1];
			$length = $stop - $start;
#			print $length."\t"; #testing
			$X12S = substr($seq, $start, $length);

			print OUT ">$id\n$X12S\n";

		}
		else {
			print "Cannot find 12S in $id: $!\n";
		}
	}
	$i++;
	$seq="";

}
$i=0;

close OUT;

#!/usr/bin/perl

#June 24, 2013 by Terri Porter
#Script to parse marker.blastn using readid_OTUcentorid.map and readid_sample.map and print custom blastn files for each marker
#usage perl parse_genbank_sort_sample.plx readid_OTUcentroid.map readid_sample.map file.blastn

use strict;
use warnings;
use Bio::SearchIO;
use Bio::SearchIO::Writer::TextResultWriter;

#declare var
my $in;
my $filename; #blastn infile
my $result;
my $queryline;
my $queryid;
my $readid;
my $OTUcentroid;
my $original;
my $new;
my $readidline;
my $j=0;
my $i=0;
my $sample;
my $outfile;
my $line;
my $sample_pooled;
my $out;
my $writer;

#declare array
my @in;
my @in2;
my @queryline;
my @line;
my @readidline;

#declare hash
my %OTUcentroid; #indexed by OTUcentroid
my %sample; #indexed by readid
my %flag; #indexed by sample_pooled

#hash OTUcentroids
open (IN, "<", $ARGV[0]) || die "Error cannot open readid_OTUcentroid.map:$!\n";
@in = <IN>;
close IN;

while($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$readid = $line[0];
	$OTUcentroid = $line[1];
	if (exists $OTUcentroid{$OTUcentroid}) {
		$original = $OTUcentroid{$OTUcentroid};
		$new = $original."|".$readid;
		$OTUcentroid{$OTUcentroid} = $new;
	}
	else {
		$OTUcentroid{$OTUcentroid} = $readid;
	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$OTUcentroid=();
}
$i=0;

#test
#while (my ($key,$value) = each (%OTUcentroid) ) {
#	print "$key => $value\n";
#}

#hash samples

open (IN2, "<", $ARGV[1]) || die "Error cannot open readid_sample.map: $!\n";
@in2 = <IN2>;
close IN2;

while ($in2[$i]) {
	$line = $in2[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$readid = $line[0];
	$sample = $line[1];
	if (exists $sample{$readid}) {
		print "Found readid $readid and sample twice\n";
	}
	else {
		$sample{$readid} = $sample;
	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$sample=();
}
$i=0;

#test
#while ( my($key,$value) = each (%sample) ) {
#	print "$key => $value\n";
#}

#parse BLAST reports
$filename = $ARGV[2];
chomp $filename;

$in = Bio::SearchIO -> new(	-file 	=> 	$filename,
							-format => 	'blast');
#print "stream:$stream\n";

while ( $result = $in -> next_result) {
	$queryline = $result->query_name;
	@queryline = split(/;/,$queryline);
	$queryid = $queryline[0];

	if (exists $OTUcentroid{$queryid}) {
		$readidline = $OTUcentroid{$queryid};
		@readidline = split(/\|/, $readidline);
		
		while ($readidline[$j]) {
			$readid = $readidline[$j];
			if (exists $sample{$readid}) {
				$sample = $sample{$readid};
				$sample =~ /((PAD|WC)\d+)\D+/;
				$sample_pooled = $1;
				if (exists $flag{$sample_pooled}) {
					$j++;
					next;
				}
				else {
					$flag{$sample_pooled} = 1;

					$outfile = $sample_pooled.".blastn";

					$writer = Bio::SearchIO::Writer::TextResultWriter-> new();
					$out = Bio::SearchIO -> new(-writer => $writer,
											-file => ">>$outfile");
					$out-> write_result($result);
				}

			}
			else {
				print "Can't find sample for $readid\n";
			}
			$j++;
		}
		$j=0;
		%flag=(); #only print OTUcentroid blast record once per sample_pooled

	}
	else {
		print "Can't find OTUcentroid that matches query $queryid\n";
	}

#	print "$queryid\n";
	
}

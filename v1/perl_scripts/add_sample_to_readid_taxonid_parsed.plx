#!/usr/bin/perl
#Oct. 10, 2013 edit to print OTUcentroids instead of each readid
#Oct. 9, 2013 by Terri Porter
#Script to add two columns for sample (PAD|WC) to target.readid.taxonid.parsed
#usage perl add_sample_to_readid_taxonid_parsed.plx cat_readid_sample.map readid_taxonid.parsed cat_readid_OTUcentroid.map

use strict;
use warnings;
use Sort::Fields;

#declare var
my $line;
my $i=0;
my $readid;
my $taxonid;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $lineage;
my $readidline;
my $centroid;
my $original;
my $new;
my $j=0;
my $pad1="";
my $pad11="";
my $pad14="";
my $pad3="";
my $pad33="";
my $pad37="";
my $pad38="";
my $pad4="";
my $wc1="";
my $wc2="";
my $wc3="";
my $wc4="";
my $wc5="";
my $wc6="";
my $wc7="";
my $wc8="";
my $sample;
my $k;
my $nextline;
my $nextreadid;
my $nextcentroid;

#declare array
my @map; #readid sample
my @in;
my @map2; #readid OTUcentroid
my @line;
my @readidline;
my @centroid;
my @nextline;
my @temp; #readid.centroid.samples.lineage
my @sorted;

#declare hash
my %map; #indexed by readid, value eq sample
my %map2; #indexed by readid value eq centroid
my %centroid; #indexed by centroid value eq readid list

open (MAP, "<", $ARGV[0]) || die "Error cannot open readid_sample.map file: $!\n";
@map = <MAP>;
close MAP;

open (IN, "<", $ARGV[1]) || die "Error cannot open target.readid.taxonid.parsed: $!\n";
@in = <IN>;
close IN;

open (CENTROID, "<", $ARGV[2]) || die "Error cannot open cat_readid_OTUcentroid.map: $!\n";
@centroid = <CENTROID>;
close CENTROID;

#hash centroid map file cat_readid_OTUcentroid.map
while ($centroid[$i]) {
	$line = $centroid[$i];
	chomp $line;

	if (length($line)>0) {
		@line = split(/\t/,$line);
		$readid = $line[0];
		$centroid = $line[1];
		if (exists $centroid{$centroid}) {
			$original = $centroid{$centroid};
			$new = $original."|".$readid;
			$centroid{$centroid} = $new;
		}
		else {
			$centroid{$centroid} = $readid;
		}
	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$centroid=();
	$original=();
	$new=();
}
$i=0;

#hash map file readid_sample.map
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	if (length($line)>0) {
		@line = split(/\t/,$line);
		$readid = $line[0];
		$sample = $line[1];
		$map{$readid}=$sample;
	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$sample=();
}
$i=0;

#parse target file
open (OUT, ">>", "centroid_samples_lineage.txt") || die "Error cannot open outfile :$!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	
	if (length($line)>0) {
		@line = split(/\t/,$line);
#		print "@line\n";
		$centroid = $line[0];
		$taxonid = $line[1];

		if ($taxonid > 0 && length($taxonid)>0 ) {
			$kingdom = $line[2];
			$phylum = $line[3];
			$class = $line[4];
			$order = $line[5];
			$family = $line[6];
			$genus = $line[7];
			$species = $line[8];
			$lineage = $kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species;

			if (exists $map{$centroid}) {
				$sample = $map{$centroid};
				$sample =~ s/^F_//;# for 16S seqs
				$sample =~ s/^R_//;
				check_sample();
			}
			else {
				print "Error cannot find sample for readid $readid\n";
			}
			
			if (exists $centroid{$centroid}) {
#				print "checking next line for centroid/readid\n";
				$readidline = $centroid{$centroid};
				@readidline = split(/\|/,$readidline);

				while ($readidline[$j]) {
					$readid = $readidline[$j];
					
					if ($map{$readid}) {
						$sample = $map{$readid};
						$sample =~ s/^F_//;
						$sample =~ s/^R_//;
						check_sample();
					}
					else {
						print "Error cannot find sample for readid $readid in readidline for centroid $centroid:$!\n";
					}
					$j++;
				}
				$j=0;
				print OUT "$centroid\t$taxonid\t$pad1\t$pad11\t$pad14\t$pad3\t$pad33\t$pad37\t$pad38\t$pad4\t$wc1\t$wc2\t$wc3\t$wc4\t$wc5\t$wc6\t$wc7\t$wc8\t$lineage\n";
			}
			else {
				print "Error cannot find readid list for centroid $centroid:$!\n";
			}
		}
	}	
	$i++;
	$line=();
	$centroid=();
	$taxonid=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
	$lineage=();
	$pad1="";
	$pad11="";
	$pad14="";
	$pad3="";
	$pad33="";
	$pad37="";
	$pad38="";
	$pad4="";
	$wc1="";
	$wc2="";
	$wc3="";
	$wc4="";
	$wc5="";
	$wc6="";
	$wc7="";
	$wc8="";
	$k=();
	$nextline=();
	@nextline=();
	$nextreadid=();
	$nextcentroid=();
}
$i=0;
close OUT;

####################

sub check_sample {

if ($sample =~ /PAD1(A|B|C)/) {
	$pad1 = 1;
}
if ($sample =~ /PAD11(A|B|C)/) {
	$pad11 = 1;
}
if ($sample =~ /PAD14(A|B|C)/) {
	$pad14 = 1;
}
if ($sample =~ /PAD3(A|B|C)/) {
	$pad3 = 1;
}
if ($sample =~ /PAD33(A|B|C)/) {
	$pad33 = 1;
}
if ($sample =~ /PAD37(A|B|C)/) {
	$pad37 = 1;
}
if ($sample =~ /PAD38(A|B|C)/) {
	$pad38 = 1;
}
if ($sample =~ /PAD4(A|B|C)/) {
	$pad4 = 1;
}
if ($sample =~ /WC1(A|B|C)/) {
	$wc1 = 1;
}
if ($sample =~ /WC2(A|B|C)/) {
	$wc2 = 1;
}
if ($sample =~ /WC3(A|B|C)/) {
	$wc3 = 1;
}
if ($sample =~ /WC4(A|B|C)/) {
	$wc4 = 1;
}
if ($sample =~ /WC5(A|B|C)/) {
	$wc5 = 1;
}
if ($sample =~ /WC6(A|B|C)/) {
	$wc6 = 1;
}
if ($sample =~ /WC7(A|B|C)/) {
	$wc7 = 1;
}
if ($sample =~ /WC8(A|B|C)/) {
	$wc8 = 1;
}

}

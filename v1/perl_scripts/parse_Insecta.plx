#!/usr/bin/perl
#May 10, 2012 by Terri Porter
#Script to parse Insecta.txt from parse_tsv.plx
#Filter by min 500 bp length and 0 ambiguities
#Create processid->BIN and fasta-formatted processid->nucraw files
#usage perl parse_Insecta.plx Insecta.txt insect_genera.txt BIN_lineage.map

use strict;
use warnings;

#declare var
my $i=0;
my $genus;
my $bin_lineage;
my $bin;
my $lineage;
my $line;
my $processid;
my $species_reg;
my $BIN;
my $nucraw;
my $length;
my $lengthCutoff=500;##### Edit this here #####
my $j=0;
my $base;
my $flag=0;
my $designation;
my $flag2=0;

#declare array
my @in;
my @in2;
my @in3;
my @bin_lineage;
my @line;
my @species_reg;
my @nucraw;
my @lineage;
my @designation;

#declare hash
my %genus;
my %bin_lineage;

open (IN_GEN, "<", $ARGV[1]) || die "Error cannot open second infile: $!\n";
@in2 = <IN_GEN>;
close IN_GEN;

#dereplicate target genera
while ($in2[$i]) {
	$genus = $in2[$i];
	chomp $genus;
	$genus{$genus} = 1;
	$i++;
}
$i=0;

open (IN_BIN, "<", $ARGV[2]) || die "Error cannot open third infile: $!\n";
@in3 = <IN_BIN>;
close IN_BIN;

while ($in3[$i]) {
	$bin_lineage = $in3[$i];
	chomp $bin_lineage;

	@bin_lineage = split(/\t/, $bin_lineage);
	$bin = $bin_lineage[0];
	$lineage = $bin_lineage[1];
	$bin_lineage{$bin} = $lineage;
	$i++;
}
$i=0;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT_BIN, ">>", "processid_BIN.map") || die "Error cannot open outfile: $!\n";

open (OUT_NUC, ">>", "processid_NUC.fasta") || die "Error cannot open outfile2: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^processid/) { #header line
		$i++;
		next;
	}
	else {
		@line = split(/\t/,$line);
		$processid = $line[0];
		#print "processid:$processid\n";#test

		$species_reg = $line[6];
		@species_reg = split(" ",$species_reg);
		$BIN = pop(@species_reg);

		print OUT_BIN "$processid\t$BIN\n";

		$nucraw = $line[7];
		#print "$nucraw\n";#test
		@nucraw = split("", $nucraw);
		$length = scalar(@nucraw);
		#print "length:$length\n";
		if ($length >= $lengthCutoff) { #lengthCutoff filter

			while ($nucraw[$j]) { #ambiguous base filter
				$base = $nucraw[$j];
				if ($base !~ /(A|C|G|T)/) {
					$flag = 1;
				}
				else {
					$j++;
					next;
				}
				$j++;
			}
			$j=0;
			#print "flag:$flag\n";#test
			if ($flag == 0) {
				$lineage = $bin_lineage{$BIN};
				#print "lineage:$lineage\n";#test
				#error during testing because haven't collected all BIN pages yet
				if ($lineage =~ /\[Genus\]/) {# need to accomodate BINS with more than one genus designation
					@lineage = split(/\|/,$lineage);
					foreach $designation (@lineage) {
						if ($designation =~ /Genus/) {
							$designation =~ s/\[Genus\]//;
							push(@designation,$designation);
						}
					}

					foreach $genus (@designation) {
						if (exists $genus{$genus}) {
							$flag2 = 1;
						}
					}

					if ($flag2 == 1) {
						print OUT_NUC ">$processid\n$nucraw\n";
					}
				}
			}
			$flag=0;
			$flag2=0;
		}
	}
	$i++;
	$line=();
	@line=();
	$processid=();
	$species_reg=();
	@species_reg=();
	$BIN=();
	$nucraw=();
	@nucraw=();
	$length=();
	$base=();
}
$i=0;
close OUT_BIN;
close OUT_NUC;

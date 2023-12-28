#!/usr/bin/perl
#June 14, 2013 by Terri Porter
#Script to turn parsed Lepidoptera parsed from .tsv file from BOLD into a fasta file like Joel's Mantodea file
#usage perl create_lep_fasta.plx iBOLD_lep.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $id;
my $bin;
my $kingdom = "Metazoa";
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $country;
my $marker;
my $seq;
my $length;
my $original;
my $lineage;

#declar array
my @in;
my @line;
my @seq;

#declare hash
my %country; #indexed by country name

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "lep.fasta") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if (length($line) > 0 ) {
#		print "found a line\n";

		@line = split(/\t/,$line);
		$id = $line[0]; #use BOLD processid instead of GB record for now
		$bin = $line[1];
		$phylum = $line[2];
		if (length($phylum) == 0) {
			$phylum = "undef_".$kingdom;
		}
		$class = $line[3];
		if (length($class) == 0) {
			$class = "undef_".$phylum;
		}
		$order = $line[4];
		if (length($order) == 0) {
			$order = "undef_".$class;
		}
		$family = $line[5];
		if (length($family) == 0) {
			$family = "undef_".$order;
		}
		$genus = $line[6];
		if (length($genus) == 0) {
			$genus = "undef_".$family;
		}
#		$species = $line[7];
		$country = $line[8];
		$marker = $line[9];
		$seq = $line[10];
#		print "seq:$seq\n";

		if ($marker =~ /COI/) {
			@seq = split(//,$seq);
			$length = scalar(@seq);
#			print "found COI marker\n";

			if ($length >= 500) { #ie. min length 500 bp
#				print "found 500 bp seq\n";

				if ($seq !~ /(B|D-F|H-S|U-Z)/i) { #ie. do not allow any missing data or non-nucleotide characters
					$lineage = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus;
					print OUT ">$id $lineage\n$seq\n";

					#keep track of country stats
					if (exists $country{$country}) {
						$original = $country{$country};
						$original++;
						$country{$country} = $original;
					}
					else {
						$country{$country} = 1;
					}

				}
				else {
					print "$id contained gaps or nucleotide ambiguities\n";
				}
			}
		}
		else {
			print "$id was a short seq\n";
		}
	}
	$i++;
	$line=();
	@line=();
	$id=();
	$bin=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
	$country=();
	$marker=();
	$seq=();
	$original=();

}
$i=0;
close OUT;

open (OUT2, ">>", "lep_country.txt") || die "Error cannot open 2nd outfile: $!\n";

while ( ($country,$original) = each (%country) ) {
	print OUT2 "$country\t$original\n";
}
close OUT2;

#!/usr/bin/perl
# Teresita M. Porter, March 14, 2023
# Turn testNBC_2.fasta into a formatted BLAST database with the optional taxonid mapping
# USAGE perl create_custom_12S_NA_Vert_BLAST_db.plx testNBC_2.fasta gb_taxid_comb.map

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $accession;
my $lineage;
my $outfasta = "12S_NA_vert.fasta";
my $seq;
my $outtxt = "12S_NA_vert.txt";
my $taxid;
my $name_txt;
my $name_class;
my $gb;
my $species;

# declare arrays
my @rdp;
my @line;
#my @acc;
my @map;
#my @map2;
#my @merged;
#my @names;

# declare hashes
my %species; # key = taxid, value = 1
#my %map2; # key = taxid, value = species
my %gb; # key = gb, value = taxid
my %acc; # key = accession, value = 1

open (RDP, "<", $ARGV[0]) || die "Error, cannot open testNBC_2.fasta: $!\n";
@rdp = <RDP>;
close RDP;

open (OUT1, ">>", $outfasta) || die "Error, cannot open outfasta: $!\n"; 

# reformat fasta file to have just the acc in the header, and uppercase seq on one line
while($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/ /, $line);
		$accession = $line[0];
		$accession =~ s/^>//g;
		$lineage = $line[1];
		$acc{$accession} = 1;
		print OUT1 ">".$accession."\n";
	}
	else {
		$seq = $line;
		$seq = uc $seq; # upper case
		print OUT1 $seq."\n";
	}

	$i++;
	$line=();
	@line=();
	$accession=();
	$lineage=();
	$seq=();

}
$i=0;
close OUT1;

# parse gb_taxid map
open (MAP, "<", $ARGV[1]) || die "Error cannot open gb_taxid_comb.map: $!\n";
@map = <MAP>;
close MAP;

while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$gb = $line[0];
	$taxid = $line[1];
	if (exists $acc{$gb}) {
		$gb{$gb} = $taxid;
	}

	$i++;
	$line=();
	@line=();
	$gb=();
	$taxid=();
}
$i=0;

open (OUT2, ">>", $outtxt) || die "Error cannot open outtxt: $!\n";

while (my ($gb, $taxid) = each %gb ) {
	print OUT2 "$gb\t$taxid\n";
}

close OUT2;


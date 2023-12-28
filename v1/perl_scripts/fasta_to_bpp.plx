#!/usr/bin/perl
#Terri Porter, August 19, 2010
#Script to parse fasta file into the phylip infile and map file needed for bpp
#Usage $perl fasta_to_table.plx genbank_fasta.txt

use strict;
use warnings;

#declare variables
my $line;
my $gi;
my $genus;
my $species;
my $details;
my $voucher;
my $i=0; #counter
my $substr;
my $l=0; #counter
#my $letter;
my $char;
my $j=0; #counter
my $k=0; #counter
my $sum;
my $average;
my $label;
my $count;
my $remaining;

#declare array
my @species;
my @chars;
my @num_char;
my @label;

#script linking subroutines
phylip_temp(); #create temporary phylip file
phylip_counts(); #get num taxa and char, limit labels to 10 char
phylip(); #print final phylip file
unlink("phylip_temp.txt");
unlink("phylip_temp2.txt");
map_temp(); #create temporary map file
map_final(); #create final map file
unlink("map_temp.txt");

#subroutine to create a temporary phylip file, add #tax and # char later
sub phylip_temp {

#open fasta file
open (IN, '<', $ARGV[0]) || die ("Error: Cannot open infile $!");

#open tempfile, phylip file without the header
open (TMP1, '>>', 'phylip_temp.txt') || die ("Error: Cannot create phylip tempfile $!");
open (TAX, '>>', 'taxon_map.txt') || die ("Error: Cannot create taxon map: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ /gi\|(\d+)\|gb\|.+\|\s{1}(\w+)\s{1}(\w+)\s{1}(.+)voucher\s{1}(.+)\s{1}/;
		$gi = $1;
		$genus = $2;
		$species = $3;
		$details = $4;
		$voucher = $5;
		$i++;
		$substr = substr $species, 0, 3;
		print TMP1 "$substr^$i\t";
		print TAX "$substr\t$i\t$gi\t$genus\t$species\t$details\t$voucher\n";
		@species=();
	}
	else {
		print TMP1 "$line\n";
	}
}
close IN;
close TMP1;
close TAX;
$i=0;
}

#subroutine to create a map_tempfile
sub map_temp {

#open fasta file
open (IN, '<', $ARGV[0]) || die ("Error: Cannot open infile $!");

#create map tempfile
open (TMP2, '>>', 'map_temp.txt') || die ("Error: Cannot create map tempfile $!");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ /gi\|(\d+)\|gb\|.+\|\s{1}(\w+)\s{1}(\w+)\s{1}(.+)voucher\s{1}(.+)\s{1}cytochrome/;
		$gi = $1;
		$genus = $2;
		$species = $3;
		$details = $4;
		$voucher = $5;
		$i++;
		print TMP2 "$i\t$species\n";
	}
	else {
		next;
	}
}
close IN;
close TMP2;
}

#subroutine to create final map file
sub map_final {

#open map_tempfile
open (IN, '<', 'map_temp.txt') || die ("Error: Cannot open map_tempfile $!");

#create map
open (MAP, '>>', 'map.txt') || die ("Error: Cannot create map file $!");

while (<IN>) {
	$line = $_;
	chomp $line;
	
	if ($line =~ s/mazaeus/A/) {
		print MAP "$line\n";
	}
	elsif ($line =~ s/polymnia/B/) {
		print MAP "$line\n";
	}
	elsif ($line =~ s/lysimnia/C/) {
		print MAP "$line\n";
	}
	elsif ($line =~ s/menapis/D/) {
		print MAP "$line\n";
	}
}
close IN;
close MAP;
}

#subroutine to count sequences and characters in phylip file
sub phylip_counts {

#open phylip_temp.txt
open (IN, '<', 'phylip_temp.txt') || die ("Error: Cannot open phylip_temp $!");
open (TMP3, '>>', 'phylip_temp2.txt') || die ("Error: Cannot open phylip_temp2: $!");

while (<IN>) {
	$line = $_;
	chomp $line;
	if ($line =~/\w+/) {
		$line =~ /(\w{3}\^\d+)\s+(\w+)/;
		$label = $1;
		$char = $2;
		@chars = split (//, $char);
		foreach (@chars) {
			$k++;
		}
		push (@num_char, $k);
		$k=0;
		$j++;
		@label = split (//,$label);
		$count = scalar(@label);
		$remaining = 10 - $count;
		print TMP3 "$label";
		while ($remaining > 0) {
			print TMP3 " ";
			$remaining--;
		}
		print TMP3 $char."\n";
	}
}

foreach (@num_char) {
	$sum += $_;
}

$average = $sum/$j;
close IN;
close TMP3;
}

#subroutine to print final phylip file
sub phylip {

#open phylip_temp.txt
open (IN, '<', 'phylip_temp2.txt') || die ("Error:Cannot open phylip_temp $!");

#create phylip file
open (PHY, '>>', 'phylip.txt') || die ("Error: cannot create final phylip file $!");

print PHY "$j\t$average\n";

while (<IN>) {
	$line = $_;
	chomp $line;
	print PHY $line."\n";
}
close IN;
close PHY;
}

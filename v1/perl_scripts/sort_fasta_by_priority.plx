#!/usr/bin/perl
# Teresita M. Porter, July 9, 2020
# Script to sort testNBC.fasta by decreasing taxonomic priority
# Linnean binomial > BOLD bin > everything else
# USAGE perl sort_fasta_by_priority.plx testNBC.fasta

use strict;
use warnings;

# vars

my $infile = $ARGV[0];
my $linnean = $infile.".linnean";
my $bold = $infile.".bold";
my $other = $infile.".other";
my $species;

my $i=0;
my $line;
my $words;
my $j;
my $seq;

# arrays
my @in;
my @line;
my @species;

# hashes

open (IN, "<", $infile) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (LIN, ">>", $linnean) || die "Error cannot open linnean outfile: $!\n";
open (BOL, ">>", $bold) || die "Error cannot open bold outfile: $!\n";
open (OTH, ">>", $other) || die "Error cannot open other outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/;/,$line);	
		$species = $line[8];
		@species = split(/_/,$species);
		$words = scalar(@species);

		if ($words == 2) { #check for Linnean binomial species name
			if ($species !~ /\d+/ && $line !~/(sp\.|cf\.|aff\.)/) { # ensure no numeric characters or species-like names
				$j = $i+1;
				$seq = $in[$j];
				chomp $seq;
				
				print LIN $line."\n".$seq."\n";
				$j=();
			}
			else {
				# contains numeric or species-like names, so not a Linnean binomial species name
				if ($species =~ /BOLD:/) { # contains a BOLD bin identifier
					$j = $i+1;
					$seq = $in[$j];
					chomp $seq;

					print BOL $line."\n".$seq."\n";
					$j=();
				}
				else {
					# lowest priority taxonomic resolution (not Linnean binomial species name, does not contain BOLD bin identifier)
					$j = $i+1;
					$seq = $in[$j];
					chomp $seq;

					print OTH $line."\n".$seq."\n";
					$j=();
				}

			}	
		}
		else {
			# more than two words so not a Linnean binomial species name

			if ($species =~ /BOLD:/) { # contains a BOLD bin identifier
				$j = $i+1;
				$seq = $in[$j];
				chomp $seq;

				print BOL $line."\n".$seq."\n";
				$j=();
			}
			else {
				# lowest priority taxonomic resolution (not Linnean binomial species name, does not contain BOLD bin identifier)
				$j = $i+1;
				$seq = $in[$j];
				chomp $seq;

				print OTH $line."\n".$seq."\n";
				$j=();
			}
		}
	}
	else {
		# skip over seq line
	}
	$i++;

}
$i=0;

close LIN;
close BOL;
close OTH;

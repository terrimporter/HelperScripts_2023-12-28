#!/usr/bin/perl
# Teresita M. Porter, March 28, 2023
# Script to filter BOLD bins
# USAGE perl filter_bins.plx bin_species.uniq > bin_species.filtered

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $i=0;
my $line;
my $bin;
my $species;
my $length;
my $first;
my $second;
my $third;

# declare arrays
my @in;
my @line;
my @species;

# declare hashes
#my %map; #key = bin, value = species

open (IN, "<", $ARGV[0]) || die "Error cannot open infile:$!\n";
@in = <IN>;
close IN;

while($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$bin = $line[0];
	$bin =~ s/ //g;
	$species = $line[1];
	$species =~ s/^\s//g; # remove white spaces
	$species =~ s/\s$//g;
#	print "$bin\t$species";

	# check for binomial species names
	@species = split(/ /, $species);
	$length = scalar @species;
	

	# skip over header
	if ($i==0) {
		$i++;
		next;
	}

	# skip over records if no BIN
	# probably not COI
	elsif ($bin =~ /^None/) {
#		print "Found None in bin field $bin\n";
		$i++;
		next;
	}
	# skip over records if no species name available
	# just keep whatever popped up from GenBank
	elsif ($species =~ /^None/) {
#		print "Found None in species field $species\n";
		$i++;
		next;
	}
	# skip over records if contains numbers
	# Linnean binomials shouldn't contain numbers 
	elsif ($species =~ /\d/) {
#		print "Found species with numbers $species\n";
		$i++;
		next;
	}
	# skip over records if uncertain ID
	elsif ($species =~ /(sp\.|sp\s|cf\.|aff\.|nr\.|\scf\s|\saff\s|\snr\s|pr\.|\spr\s|prob\.|cfr\.|gr\.|n\.sp|\/)/) {
#		print "Found species that are not sufficiently identified $species\n";
		$i++;
		next;
	}
	# not a binomial species name
	elsif ($length < 2) {
		$i++;
		next;
	}
	# handle hybrids, just keep them
	elsif ($species =~ / x /i) {
		$species =~ s/ /_/g;
		print "$bin\t$species\n";
	}
	# try to salvage a binomial species name
	elsif ($length > 2) {
		$first = $species[0];
		$second = $species[1];
		$third = $species[2];
		if ($second =~ /(\(|\))/) {
			$second = $third;
			$species = $first."_".$second;
			print "$bin\t$species\n";
#			print "found punctuation in first species word $species\n";
		}
#		print "found species with more than 2 fields $species\n";
	}
	else {
#	$map[$bin] = $species;
		$species =~ s/ /_/g;
		print "$bin\t$species\n";
	}

	$i++;
	$line=();
	$bin=();
	$species=();

}
$i=0;

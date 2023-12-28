#!/usr/bin/perl
#May 25, 2012 by Terri Porter
#Script to parse lichen.xls into fields (synonyms and non-lichens removed by hand)
#usage perl parse_lichen.plx lichen.xls

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $genus;
my $Genus;
my $species;
my $key;
my $value;
my $name;

#declare array
my @in;
my @line;
my @line2;
my @name;

#declare hash
my %name;

open (IN,"<",$ARGV[0]) || die "Error cannot open lichen.xls";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^[A-Z]/) { #check for UPPERCASE genus name
		@line = split(" ", $line);
		$genus = $line[0];
		#print "$genus\t";
		$genus = lc($genus); #change to lowercase
		$Genus = ucfirst($genus); #capitalize first letter only
		#print "$Genus\n";
		@line=();
	}
	elsif ($line =~ /^[a-z]/) { #check for lowercase species name
		@line2 = split(" ", $line);
		$species = $line2[0];
		$key = $Genus."|".$species;
		$name{$key} = 1;
		@line2=();
	}

	$i++;
	$line=();
}
$i=0;

open (OUT, ">>", "lichen.txt") || die "Error cannot open lichen.txt: $!\n";

while (($key,$value) = each (%name)) {
	$name = $key;
	@name = split(/\|/,$name);
	$genus = $name[0];
	$species = $name[1];
	print OUT "$genus\t$species\n";

	$name=();
	@name=();
	$genus=();
	$species=();
}
close OUT;

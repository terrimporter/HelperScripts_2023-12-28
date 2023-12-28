#!/usr/bin/perl
#Oct.6, 2011 by Terri Porter
#Script to merge 'best' files from parse_clone_html_files_rank4.plx scripts
#usage perl merge_sap_best_results.plx clone_summary_best_kingdom.txt clone_summary_best_phylum.txt
#clone_summary_best_class.txt clone_summary_best_order.txt clone_summary_best_family.txt
#clone_summary_best_genus.txt clone_summary_best_species.txt gi.query

use strict;
use warnings;

#declare var
my $cutoff;
my $a=0;
my $line;
my $gi;
my $kingdom;
my $stat;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;

#declare array
my @line;
my @kingdom;
my @phylum;
my @class;
my @order;
my @family;
my @genus;
my @species;
my @gi;

#declare hash
my %kingdom;
my %phylum;
my %class;
my %order;
my %family;
my %genus;
my %species;

#####edit value used for statistical support cutoff here#####
$cutoff = 0;

open (KINGDOM,"<",$ARGV[0]) || die "Error cannot read from clone_summary_best_kingdom.txt: $!\n";
@kingdom = <KINGDOM>;
close KINGDOM;

while ($kingdom[$a]) {
	$line = $kingdom[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$kingdom = $line[1];
		$stat = $line[2];
		if ($stat >= $cutoff ) {
			$kingdom{$gi} = $kingdom;
		}
	}
	$a++;
	@line=();
	$gi=();
	$kingdom=();
	$stat=();
}
$a=0;

open (PHYLUM,"<",$ARGV[1]) || die "Error cannot read from clone_summary_best_phylum.txt: $!\n";
@phylum = <PHYLUM>;
close PHYLUM;

while ($phylum[$a]) {
	$line = $phylum[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$phylum = $line[1];
		$stat = $line[2];
		if ($stat >= $cutoff ) {
			$phylum{$gi} = $phylum;
		}
	}
	$a++;
	@line=();
	$gi=();
	$phylum=();
	$stat=();
}
$a=0;

open (CLASS,"<",$ARGV[2]) || die "Error cannot read from clone_summary_best_class.txt: $!\n";
@class = <CLASS>;
close CLASS;

while ($class[$a]) {
	$line = $class[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$class = $line[1];
		$stat = $line[2];
		if ($stat >= $cutoff ) {
			$class{$gi} = $class;
		}
	}
	$a++;
	@line=();
	$gi=();
	$class=();
	$stat=();
}
$a=0;

open (ORDER,"<",$ARGV[3]) || die "Error cannot read from clone_summary_best_order.txt: $!\n";
@order = <ORDER>;
close ORDER;

while ($order[$a]) {
	$line = $order[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$order = $line[1];
		$stat = $line[2];
		if ($stat >= $cutoff ) {
			$order{$gi} = $order;
		}
	}
	$a++;
	@line=();
	$gi=();
	$order=();
	$stat=();
}
$a=0;

open (FAMILY,"<",$ARGV[4]) || die "Error cannot read from clone_summary_best_family.txt: $!\n";
@family = <FAMILY>;
close FAMILY;

while ($family[$a]) {
	$line = $family[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$family = $line[1];
		$stat = $line[2];
		if ($stat >= $cutoff ) {
			$family{$gi} = $family;
		}
	}
	$a++;
	@line=();
	$gi=();
	$family=();
	$stat=();
}
$a=0;

open (GENUS,"<",$ARGV[5]) || die "Error cannot read from clone_summary_best_genus.txt: $!\n";
@genus = <GENUS>;
close GENUS;

while ($genus[$a]) {
	$line = $genus[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$genus = $line[1];
		$stat = $line[2];
		if ($stat >= $cutoff ) {
			$genus{$gi} = $genus;
		}
	}
	$a++;
	@line=();
	$gi=();
	$genus=();
	$stat=();
}
$a=0;

open (SPECIES,"<",$ARGV[6]) || die "Error cannot read from clone_summary_best_species.txt: $!\n";
@species = <SPECIES>;
close SPECIES;

while ($species[$a]) {
	$line = $species[$a];
	chomp $line;
	if ($line !~ /nil\tnil/) {
		@line = split(/\t/,$line);
		$gi = $line[0];
		$species = $line[2];
		$stat = $line[3];
		if ($stat >= $cutoff ) {
			if ($species !~ /(sp\.|cf\.|aff\.)/) {
				$species{$gi} = $species;
			}
		}
	}
	$a++;
	@line=();
	$gi=();
	$species=();
	$stat=();
}
$a=0;

open (GI,"<",$ARGV[7]) || die "Error cannot read from gi.query: $!\n";
@gi = <GI>;
close GI;

open (OUT,">>","merged_sap.txt") || die "Error cannot write to merged_sap.txt: $!\n";

while ($gi[$a]) {
	$gi = $gi[$a];
	chomp $gi;

	$kingdom = $kingdom{$gi};
	$phylum = $phylum{$gi};
	$class = $class{$gi};
	$order = $order{$gi};
	$family = $family{$gi};
	$genus = $genus{$gi};
	$species = $species{$gi};

	if ($kingdom !~ /^\S+/) {
		$kingdom = "nil";
	}
	if ($phylum !~ /^\S+/) {
		$phylum = "nil";
	}
	if ($class !~ /^\S+/) {
		$class = "nil";
	}
	if ($order !~ /^\S+/) {
		$order = "nil";
	}
	if ($family !~ /^\S+/) {
		$family = "nil";
	}
	if ($genus !~ /^\S+/) {
		$genus = "nil";
	}
	if ($species !~ /^\S+/) {
		$species = "nil";
	}
	
	if (($kingdom and $phylum and $class and $order and $family and $genus and $species) eq "nil") {
		$a++;
		next;
	}
	else {
		print OUT "$gi\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
	}
	$a++;
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
}
$a=0;

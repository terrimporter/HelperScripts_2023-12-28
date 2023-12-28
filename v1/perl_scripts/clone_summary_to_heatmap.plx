#!/usr/bin/perl
#Oct.6,2010 by Terri Porter
#Script to convert clone_summary.mapped to a tabular form suitable for making a heatmap
#usage $perl clone_summary_to_heatmap.plx clone_summary.mapped

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $x;
my $i=0;
my $line;
my $speciesQ;
my $speciesA;
my $pp;
my $newline;
my $means;
my $mean_bayanus;
my $mean_cariocanus;
my $mean_cerevisiae;
my $mean_kudriavzevii;
my $mean_mikatae;
my $mean_paradoxus;
my $mean_pastorianus;

#declare array
my @clone_summary;
my @line;
my @cerevisiae;
my @bayanus;
my @cariocanus;
my @kudriavzevii;
my @mikatae;
my @paradoxus;
my @pastorianus;
my @bayanus_pp;
my @cariocanus_pp;
my @cerevisiae_pp;
my @kudriavzevii_pp;
my @mikatae_pp;
my @paradoxus_pp;
my @pastorianus_pp;

open (IN,"<",$ARGV[0]) || die ("Error:$!\n");
@clone_summary=<IN>;
close IN;

foreach $x (@clone_summary) {
	chomp $x;
}

while ($clone_summary[$i]) {
	$line = $clone_summary[$i];
	@line = split(/\t/,$line);
	$speciesQ = $line[2];
	$speciesA = $line[4];
	$pp = $line[5];
	$newline = "$speciesA\t$pp";
	if ($speciesQ =~ /cerevisiae/) {
		push(@cerevisiae,$newline);
	}
	elsif ($speciesQ =~ /bayanus/) {
		push(@bayanus,$newline);
	}
	elsif ($speciesQ =~ /cariocanus/) {
		push (@cariocanus,$newline);
	}
	elsif ($speciesQ =~ /kudriavzevii/) {
		push (@kudriavzevii,$newline);
	}
	elsif ($speciesQ =~ /mikatae/) {
		push (@mikatae,$newline);
	}
	elsif ($speciesQ =~ /paradoxus/) {
		push (@paradoxus,$newline);
	}
	elsif ($speciesQ =~ /pastorianus/) {
		push (@pastorianus,$newline);
	}
	$i++;
}
$i=0;
empty_arrays();
while ($bayanus[$i]) {
	$line = $bayanus[$i];
	@line = split(/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];	
	compile_pps();
	$i++;
}
$i=0;
get_averages();
print "\tbayanusA\tcariocanusA\tcerevisiaeA\tkudriavzeviiA\tmikataeA\tparadoxusA\tpastorianusA\n";
print "bayanusQ\t$means\n";
empty_arrays();
while ($cariocanus[$i]){
	$line = $cariocanus[$i];
	@line = split (/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];
	compile_pps();
	$i++;
}
$i=0;
get_averages();
print "cariocanusQ\t$means\n";
empty_arrays();
while ($cerevisiae[$i]) {
	$line = $cerevisiae[$i];
	@line = split (/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];
	compile_pps();
	$i++;
}
$i=0;
get_averages();
print "cerevisiaeQ\t$means\n";
empty_arrays();
while ($kudriavzevii[$i]) {
	$line = $kudriavzevii[$i];
	@line = split (/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];
	compile_pps();
	$i++;
}
$i=0;
get_averages();
print "kudriavzeviiQ\t$means\n";
empty_arrays();
while ($mikatae[$i]) {
	$line = $mikatae[$i];
	@line = split(/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];
	compile_pps();
	$i++;
}
$i=0;
get_averages();
print "mikataeQ\t$means\n";
empty_arrays();
while ($paradoxus[$i]){
	$line = $paradoxus[$i];
	@line = split(/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];
	compile_pps();
	$i++;
}
$i=0;
get_averages();
print "paradoxusQ\t$means\n";
empty_arrays();
while ($pastorianus[$i]){
	$line = $pastorianus[$i];
	@line = split(/\t/,$line);
	$speciesA = $line[0];
	$pp = $line[1];
	compile_pps();
	$i++;
}
get_averages();
print "pastorianusQ\t$means\n";


###################

sub empty_arrays {

	@bayanus_pp=();
	@cariocanus_pp=();
	@cerevisiae_pp=();
	@kudriavzevii_pp=();
	@mikatae_pp=();
	@paradoxus_pp=();
	@pastorianus_pp=();

}

####################

sub compile_pps {

	if ($speciesA =~ /bayanus/) {
		push(@bayanus_pp,$pp);
	}
	if ($speciesA =~ /cariocanus/) {
		push(@cariocanus_pp,$pp);
	}
	if ($speciesA =~ /cerevisiae/) {
		push(@cerevisiae_pp,$pp);
	}
	if ($speciesA =~ /kudriavzevii/) {
		push (@kudriavzevii_pp,$pp);
	}
	if ($speciesA =~ /mikatae/) {
		push(@mikatae_pp,$pp);
	}
	if ($speciesA =~ /paradoxus/) {
		push(@paradoxus_pp,$pp);
	}
	if ($speciesA =~ /pastorianus/) {
		push(@pastorianus_pp,$pp)
	}

}

#####################

sub get_averages {

	my $j=0;
	$mean_bayanus = mean (@bayanus_pp);
	if (!$mean_bayanus) {
		$mean_bayanus = 0;
	}
	$mean_cariocanus = mean (@cariocanus_pp);
	if (!$mean_cariocanus) {
		$mean_cariocanus = 0;
	}
	$mean_cerevisiae = mean (@cerevisiae_pp);
	if (!$mean_cerevisiae) {
		$mean_cerevisiae = 0;
	}
	$mean_kudriavzevii = mean (@kudriavzevii_pp);
	if (!$mean_kudriavzevii) {
		$mean_kudriavzevii = 0;
	}
	$mean_mikatae = mean (@mikatae_pp);
	if (!$mean_mikatae) {
		$mean_mikatae = 0;
	}
	$mean_paradoxus = mean(@paradoxus_pp);
	if (!$mean_paradoxus) {
		$mean_paradoxus = 0;
	}
	$mean_pastorianus = mean(@pastorianus_pp);
	if (!$mean_pastorianus) {
		$mean_pastorianus = 0;
	}
	$means = "$mean_bayanus\t$mean_cariocanus\t$mean_cerevisiae\t$mean_kudriavzevii\t$mean_mikatae\t$mean_paradoxus\t$mean_pastorianus";

}
	

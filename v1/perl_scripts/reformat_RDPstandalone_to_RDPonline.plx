#!/usr/bin/perl

#Aug.12,2016 by Terri Porter.  Script to reformat RDP standalone output into a format that can be imported by MEGAN as if it were from from RDP online classifier.
#USAGE perl reformat_RDPstandalone_to_RDPonline.plx DT1.fasta.out

use strict;
use warnings;

#var
my $infile;
my $outfile;
my $i=0;
my $line;
my $field4;
my $field7;
my $field10;
my $field13;
my $field16;
my $field19;
my $newid;
my $strand;
my $kingdom;
my $kingdom_bs;
my $phylum;
my $phylum_bs;
my $class;
my $class_bs;
my $order;
my $order_bs;
my $family;
my $family_bs;
my $genus;
my $genus_bs;
my $newline;

#array
my @NBC;
my @split;
my @newline;

open (NBC, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@NBC = <NBC>;
close NBC;

$infile = $ARGV[0];
$infile =~ s/\.out//g;
$outfile = $infile."_classified.txt";

open (OUT, ">>", $outfile) || die "Cannot open outfile: $!\n";

print OUT "Classifier: RDP Naive Bayesian rRNA Classifier Version 2.12.\nTaxonomical Hierarchy: RDP COIv2 training set\nQuery File: XXX.fasta\nSubmit Date: Fri Aug 12 09:52:36 EDT 2016\nConfidence threshold (for classification to Root ONLY): 90%\nSymbol +/- indicates predicted sequence orientation\n\n";

while ($NBC[$i]) { #parse each OTU
	$line = $NBC[$i];
	chomp $line;

	@split = split ("\t", $line);

	$newid = $split[0];
	$newid =~ s/;$//g;
	$newid =~ s/;/_/g;
	$newid =~ s/=/_/g;
	$newid =~ s/-/_/g;
		
	$strand = "+";

	$kingdom = $split[2];
		
	$field4 = $split[4];
	$kingdom_bs = $field4*100;

	$phylum = $split[5];
		
	$field7 = $split[7];
	$phylum_bs = $field7*100;

	$class = $split[8];

	$field10 = $split[10];
	$class_bs = $field10*100;

	$order = $split[11];

	$field13 = $split[13];
	$order_bs = $field13*100;

	$family = $split[14];

	$field16 = $split[16];
	$family_bs = $field16*100;

	$genus = $split[17];

	$field19 = $split[19];
	$genus_bs = $field19*100;

	$newline = $newid.";".$strand.";".$kingdom.";".$kingdom_bs."%;".$phylum.";".$phylum_bs."%;".$class.";".$class_bs."%;".$order.";".$order_bs."%;".$family.";".$family_bs."%;".$genus.";".$genus_bs."%";

	print OUT $newline."\n";

	@split=();
	$field4 = ();
	$field7 = ();
	$field10 = ();
	$field13 = ();
	$field16 = ();
	$field19 = ();
	$strand = ();
	$kingdom = ();
	$kingdom_bs = ();
	$phylum = ();
	$phylum_bs = ();
	$class = ();
	$class_bs = ();
	$order = ();
	$order_bs = ();
	$family = ();
	$family_bs = ();
	$genus = ();
	$genus_bs = ();
	$newline = ();
	$i++;
	
}
$i=0;

close OUT;

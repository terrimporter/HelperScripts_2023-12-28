#!/usr/bin/perl
#April 29, 2013 edited to grab order rank stats
#March 25, 2013 by Terri Porter
#Script to check classification tests using subusets of testNBC.fasta classified using testNBC.taxonomy
#usage perl check_NBC.plx testquery.out testNBC.fasta
#add a filter to only look at assignments with > x% bootstrap confidence

use strict;
use warnings;

#declare vars
my $i=0;
my $line;
my $scalar;
my $genusField;
my $id;
my $asst;
my $idstring;
my $correct=0; 
my $max;
my $known;
my $wrong=0;
my $comment;
my $bootstrapCutoff = 0.9;
my $bootstrapField;
my $confidence;

#declare arrays
my @nbc;
my @line;
my @known;
my @wrong;
my @notClassified;

#declare hashes
my %asst;
my %known;
my %confidence;

open (NBC, "<", $ARGV[0]) || die "Error cannot open testquery.out: $!\n";
@nbc = <NBC>;
close NBC;

open (KNOWN, "<", $ARGV[1]) || die "Error cannot open testNBC.fasta: $!\n";
@known = <KNOWN>;
close KNOWN;

#hash assts
while ($nbc[$i]) {
	$line = $nbc[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$scalar = scalar(@line);
	$genusField = $scalar-9; #edit to grab order asst
	$bootstrapField = $scalar-7; #edit to grab order bs
	$id = $line[0];
	$asst = $line[$genusField];
	$confidence = $line[$bootstrapField];

	$asst{$id} = $asst;
	$confidence{$id} = $confidence;

	$i++;
	$line=();
	@line=();
	$scalar=();
	$genusField=();
	$bootstrapField=();
	$confidence=();
	$id=();
	$asst=();
}
$i=0;

#hash known
while ($known[$i]) {
	$line = $known[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(";", $line);
		$idstring = $line[0];
		$idstring =~ /^>(\S+)\s+/;
		$id = $1;
		$scalar = scalar(@line);
		$max = $scalar-3; #edit to grab order rank
		$known = $line[$max];
		
		$known{$id} = $known;
	}
	$i++;
	$line=();
	@line=();
	$idstring=();
	$id=();
	$scalar=();
	$max=();
	$known=();
}
$i=0;

#compare asst with known
while ( ($id, $asst) = each %asst) {

	if (exists $confidence{$id}) {
		$confidence = $confidence{$id};
		
		if ($confidence > $bootstrapCutoff) {

			if (exists $known{$id}) {
				$known = $known{$id};

				if ($asst eq $known) {
					$correct++;
				}
				else {
					$wrong++;
					if (exists $confidence{$id}) {
						$confidence = $confidence{$id};
						$comment = $id."|".$asst."|".$confidence."|".$known;
						push(@wrong, $comment);
					}
					else {
						print "Cannot find confidence value for $id\n";
						$comment = $id."|".$asst."|"."|".$known;
						push(@wrong, $comment);
					}
				}
			}
			else {
				print "Error couldn't find $id among known references\n";
			}
		}
		else {
			push(@notClassified, $id);	
		}
	}
	else {
		print "Error couldn't find $id among known reference bootstraps\n";
	}
	$known=();
	$confidence=();
	$comment=();
}

$scalar = scalar(@notClassified);

print "\nUsing a taxonomic bootstrap confidence cutoff of $bootstrapCutoff:\n";
print "$correct correct assignments\n$wrong incorrect assignments\n$scalar not classified\n\n";

while ($wrong[$i]) {
	$line = $wrong[$i];
	print "WRONG:$line\n";
	$i++;
}
$i=0;

while ($notClassified[$i]) {
	$line = $notClassified[$i];
	print "NOTCLASSIFIED:$line\n";
	$i++;
}
$i=0;
print "\n";

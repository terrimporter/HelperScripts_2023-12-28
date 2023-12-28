#!/usr/bin/perl
#May 5, 2013 edit to use variable bootstrap support cutoffs for different ranks
#May 3, 2013 by Terri Porter
#Script to compare malaise morphological identifications (order) with nbc v2 and v3 classifications
#usage perl compare_morphid.plx 1052targets.fas malaise_v2.out

use warnings;
use strict;

#declare var
my $i=0;
my $line;
my $id; #morphids
my $orderMorph;
my $order; #nbc assigned
my $family;
my $genus;
my $orderVal;
my $familyVal;
my $genusVal;
my $string;
my $match=0;
my $mismatch=0;
my $genus_improved=0;
my $family_improved=0;
my $recheck=0;
my $notclassified=0;

#declare array
my @fas;
my @nbc;
my @line;
my @string;

#declare hash
my %morph; #indexed by morphid
my %asst; #indexed by morphid

open (FAS, "<", $ARGV[0]) || die "Error cannot open fasta: $!\n";
@fas = <FAS>;
close FAS;

open (NBC, "<", $ARGV[1]) || die "Error cannot open malaise.out: $!\n";
@nbc = <NBC>;
close NBC;

#hash morphids
while ($fas[$i]) {
	$line = $fas[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$id = $1;
#		print "id:$id\n";
		$id =~ /(\w+)JGP/;
		$orderMorph = $1;
#		print "orderMorph:$orderMorph\n";
		$morph{$id} = $orderMorph;
	}
	$i++;
	$line=();
	$id=();
	$orderMorph=();
}
$i=0;

#hash nbc assignments
while($nbc[$i]) {
	$line = $nbc[$i];
	chomp $line;

	@line = split(' ',$line);
	$id = $line[0];
	$order = $line[10];
	$orderVal = $line[12];
	$family = $line[13];
	$familyVal = $line[15];
	$genus = $line[16];
	$genusVal = $line[18];
	$string = $order."|".$orderVal."|".$family."|".$familyVal."|".$genus."|".$genusVal;
#	print "string:$string\n";
	$asst{$id} = $string;
	
	$i++;
	$line=();
	@line=();
	$id=();
	$order=();
	$orderVal=();
	$family=();
	$familyVal=();
	$genus=();
	$genusVal=();
	$string=();
}
$i=0;

#check if morphid order matches assigned order
#if matches, see if nbc could improve the assignment (90% or greater only)
open (OUT, ">>", "matches.txt") || die "Error couldn't open outfile matches: $!\n";
open (OUT2, ">>", "mismatches.txt") || die "Error couldn't open outfile mismatches: $!\n";
open (OUT3, ">>", "genus_improved.txt") || die "Error couldn't open outfile genus improved: $!\n";
open (OUT4, ">>", "family_improved.txt") || die "Error couldn't open outfile family improved: $!\n";
open (OUT5, ">>", "morphids_to_recheck.txt") || die "Error couldn't open outfile rechecks: $!\n";

while( ($id,$orderMorph) = each(%morph) ) {

	if (exists $asst{$id}) {
		$string = $asst{$id};
		@string = split(/\|/,$string);
		$order = $string[0];
		$orderVal = $string[1];
		$family = $string[2];
		$familyVal = $string[3];
		$genus = $string[4];
		$genusVal = $string[5];

		if ($orderMorph eq $order && $orderVal >= 0.20) {
				$match++;
				print OUT $id."|".$orderMorph."|".$string."\n";

			if ($genusVal >= 0.70) {
				$genus_improved++;
				print OUT3 $id."|".$orderMorph."|".$string."\n";
			}

			if ($familyVal >= 0.40) {
				$family_improved++;
				print OUT4 $id."|".$orderMorph."|".$string."\n";
			}
		}
		elsif ($orderMorph ne $order && $orderVal >= 0.20) {
			$mismatch++;
			print OUT2 $id,"|",$orderMorph."|".$string."\n";
#			if ($orderVal >= 0.20) {
#				$recheck++;
#				print OUT5 $id."|".$orderMorph."|".$string."\n";
#			}

		}
		else {
			$notclassified++;
		}

	}
	else {
		print "Could not find an assignment for $id\n";
	}
}

close OUT;
close OUT2;
close OUT3;
close OUT4;
close OUT5;

print "There were $match matches between the morphid and assigned order (20%+)\n";
print "There were $mismatch mismatches between the morphid and assigned order (20%+)\n";
print "$genus_improved assignments were improved to the genus rank with nbc (70%+)\n";
print "$family_improved assignments were improved to the family rank with nbc (40%+)\n";
#print "Suggest that $recheck morphids be checked\n";
print "There were $notclassified queries not classified with any significant bootstrap support\n";

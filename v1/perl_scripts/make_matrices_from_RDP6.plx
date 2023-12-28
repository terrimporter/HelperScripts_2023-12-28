#!/usr/bin/perl
#Dec. 19, 2016 edit to include all barcoding genera (from INV samples) not just morph genera (and family_undef)
#Dec.12, 2016 by Terri Porter
#Script to make presence-absence and numOccurrence matrices for each sample from RDP classifier CO1v2 output
#Customize to process Blackbrook headers, genus assignments bp>0.6, else if family bp > 0.3 then family_undef
#USAGE perl make_matrices_from_RDP6.plx 2014_mapping.txt YEAR_rdp.out

use strict;
use warnings;

#declare var
my $sampleLine;
my $sampleName;
my $cutoff;
my $minsizecutoff = 3; #use 3 to exclude singletons and doubletons
my $genusMorph;
my $line;
my $genus;
my $header;
my $site; #field site replicate name
my $size; #OTUsize, i.e. number of reads clustered in OTU
my $genusBP; 
my $genusBPcutoff = 0.6; #USE 0.6 with CO1v2 RDP training set (BE and F230 ~ 200bp )
my $family;
my $familyBP; 
my $familyBPcutoff = 0.3; #USE 0.3 cutoff with CO1v2 RDP training set then make genus asst family_undef
my $oldvalue;
my $newvalue;
my $oldvalue2;
my $newvalue2;
my $counter=0;
my $newcounter;
my $percentOccurrence;
my $i=0;
my $genusNEW; #for second round of parsing so that original genus can be retained for checking against genusTracker
my $taxon; #famil_undef

#declare array
my @sampleMappings;
my @sampleLine;
my @parsed;
my @line;
my @header;
my @sites;

#declare hash
my %sampleMappings; #key = SampleName, value = PooledRepOriginalSampleName
my %samples; #key = SampleName, value = 0 to increment for totalOccurrence per sample
my %genusTracker; #key = genus (with a BP > 0.6), value = "0"; on second round of parsing, use to check to see if genus exists, in which case do not create family_undef
my %taxonTracker; #key = genus (genusBP>0.6) or family_undef (genusBP<0.6 && familyBP>0.3, only if genus doesn't already exist)

open (IN, "<", $ARGV[0]) || die "Cannot open infile2: $!\n";
@sampleMappings = <IN>;
close IN;

foreach $sampleLine (@sampleMappings) {
	chomp $sampleLine;
	@sampleLine = split(/\t/,$sampleLine);
	$sampleName = $sampleLine[0];
	$samples{$sampleName} = "0";
}

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@parsed = <IN2>;
close IN2;

open (OUT, ">>", "totalOccurrencePerSample.txt") || die "Cannot open outfile1: $!\n";
print OUT "BarcodingTaxa";

foreach $sampleName (sort(keys(%samples))) {
	chomp $sampleName;
	print OUT "\t".$sampleName;
}
print OUT "\n";

open (OUT2, ">>", "presenceAbsencePerSample.txt") || die "Cannot open outfile2: $!\n";
print OUT2 "BarcodingTaxa";

foreach $sampleName (sort(keys(%samples))) {
	chomp $sampleName;
	print OUT2 "\t".$sampleName;
}
print OUT2 "\n";

#FIRST PASS ONLY FOR GENUSBP >= 0.6
while ($parsed[$i]) {
	$line = $parsed[$i];
	chomp $line;
	@line = split(/\t/,$line);

	$sampleLine = $line[0];
	@sampleLine = split(/;/,$sampleLine);
	$site = $sampleLine[1];
	
	if (exists $samples{$site}) {  ### only record genera from INV samples

		$genus = $line[14]; 
		$genusBP = $line[15]; 

		#check for minimum bootstrap support cutoff value, FIRST PASS, only process if genusBP >= 0.6
		if ($genusBP >= $genusBPcutoff) {
			$genusTracker{$genus} = "0"; #use to check against during SECOND round of parsing
			$taxonTracker{$genus} = "0"; #use to check against during THIRD round of parsing
		}
		else {
			$i++;
			next;
		}
	}
	else {
		$i++;
		next;
	}

	$i++;
}
$i=0;

#SECOND PASS ONLY FOR GENUSBP < 0.6 && FAMILYBP >= 0.3 to add family_undef to %taxonTracker 
while ($parsed[$i]) {
	$line = $parsed[$i];
	chomp $line;
	@line = split(/\t/,$line);

	$sampleLine = $line[0];
	@sampleLine = split(/;/,$sampleLine);
	$site = $sampleLine[1];

	if (exists $samples{$site}) { ### only record genera & family_undef from INV samples

		$family = $line[12];
		$familyBP = $line[13];
		$genus = $line[14]; 
		$genusBP = $line[15]; 

		#only process if genusBP < 0.6 && familyBP > 0.3 && genusBP doesn't already exist
		if ($genusBP < $genusBPcutoff && $familyBP >= $familyBPcutoff) {
			if (exists $genusTracker{$genus}) {
				$i++;
				next;
			}
			else {
				$taxon = $family."_undef";
				$taxonTracker{$taxon} = "0"; #use to check against during THIRD round of parsing
			}
		}
		else {
			$i++;
			next;
		}
	}
	else {
		$i++;
		next;
	}

	$i++;
}
$i=0;


### THIRD PASS ###
foreach $taxon (sort keys %taxonTracker) { ### considers all genera and family_undef not just morphGenera
	chomp $taxon;

	while ($parsed[$i]) {
		$line = $parsed[$i];
		chomp $line;
		@line = split(/\t/,$line);

		$sampleLine = $line[0];
		@sampleLine = split(/;/,$sampleLine);
		$site = $sampleLine[1];
		$size = $sampleLine[2];

		if (exists $samples{$site}) { ### only record genera & family_undef from INV samples

			$genus = $line[14];
			$genusBP = $line[15];

			if ($genusBP < $genusBPcutoff) {
				$family = $line[12]; 	
				$familyBP = $line[13];
		
				if ($familyBP >= $familyBPcutoff) {
					$genusNEW = $family."_undef";
				}
				else {
					$genusNEW = "undef"; 
				}
			}

#			$header = $line[0];
#			@header = split(/;/,$header);
#			$site = $header[1];
#			$size = $header[2];
		
#		if (exists $samples{$site}) {
		
			if ($genus eq $taxon || $genusNEW eq $taxon) {

				if ($genusBP >= $genusBPcutoff) {
					
					if ($size >= $minsizecutoff) {
			
						$oldvalue = $samples{$site}; #update samples counter for total Occurrence
						$newvalue = $oldvalue + 1;
						$samples{$site} = $newvalue;

					}
					else {
						$i++;
						next;
					}
				}
				elsif ($genusBP < $genusBPcutoff) {
					if ($familyBP >= $familyBPcutoff) {
						if (exists $taxonTracker{$genus}) {
							$i++;
							next;
						}
						else {
							if ($size >= $minsizecutoff) {
			
								$oldvalue = $samples{$site}; #update samples counter for total Occurrence
								$newvalue = $oldvalue + 1;
								$samples{$site} = $newvalue;

							}
							else {
								$i++;
								next;
							}
							
						}
					}
					else {
						$i++;
						next;
					}
				}
			}
			else {
				$i++;
				next;
			}
		}
		else {
			$i++;
			next;
		}

		$i++;
	}
		
	#totalOccurrencePerSample
	print OUT "$taxon";
	#print values for each rep
	foreach $sampleName (sort (keys (%samples))) {
		$newvalue = $samples{$sampleName};
		print OUT "\t".$newvalue;
	}
	print OUT "\n";

	#presenceAbsencePerSample
	print OUT2 "$taxon";
	foreach $sampleName (sort (keys (%samples))) {
		if ($samples{$sampleName} >= 1) {
			print OUT2 "\t1";
		}
		else {
			print OUT2 "\t0";
		}
		$samples{$sampleName} = "0"; #reset
	}
	print OUT2 "\n";

	$i=0;
}


close OUT;
close OUT2;

#!/usr/bin/perl
#Dec. 20, 2016 modified to track all barcoding genera (from INV samples) not just morphGenera
#Dec.9, 2016 by Terri Porter
#Script to make presence-absence and numOccurrence matrices for each sample and for pooled samples from RDP classifier CO1v2 output
#Customize to process Blackbrook headers, genus assignments bp>0.6, else if family bp > 0.3 then family_undef
#USAGE perl make_matrices_from_RDP3.plx sampleMapping.txt YEAR_rdp.out

use strict;
use warnings;

#declare var
my $sampleLine;
my $sampleName;
my $samplePooledName;
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
my $numReps = 3; #already known for Blackbrook 2015
my $i=0;
my $genusNEW; #for second round of parsing so that original genus can be retained for checking against genusTracker
my $taxon;

#declare array
my @genusMorph;
my @sampleMappings;
my @sampleLine;
my @parsed;
my @line;
my @header;
my @sites;

#declare hash
my %sampleMappings; #key = SampleName, value = PooledRepOriginalSampleName
my %samples; #key = SampleName, value = 0 to increment for totalOccurrence per sample
my %pooledSamples; #key = PooledRepOriginalSampleName, value = 0 to increment for totalOccurrence per samplePooled
my %genusTracker; #key = genus (with a BP > 0.6), value = "0"; on second round of parsing, use to check to see if genus exists, in which case do not create family_undef
my %taxonTracker; #key = genus (genusBP>=0.6) or family_undef (genusBP<0.6 && familyBP>0.3, only if genus doesn't already exist)

#open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
#@genusMorph = <IN>;
#close IN;

open (IN, "<", $ARGV[0]) || die "Cannot open infile2: $!\n";
@sampleMappings = <IN>;
close IN;

foreach $sampleLine (@sampleMappings) {
	chomp $sampleLine;
	@sampleLine = split(/\t/,$sampleLine);
	$sampleName = $sampleLine[0];
	$samplePooledName = $sampleLine[1];
	$sampleMappings{$sampleName} = $samplePooledName;
	$samples{$sampleName} = "0";
	$pooledSamples{$samplePooledName} = "0";
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

open (OUT3, ">>", "totalOccurrencePerPooledSample.txt") || die "Cannot open outfile3: $!\n";
print OUT3 "BarcodingTaxa";

foreach $samplePooledName (sort (keys (%pooledSamples))) {
	print OUT3 "\t".$samplePooledName;
}
print OUT3 "\n";

open (OUT4, ">>", "presenceAbsencePerPooledSample.txt" ) || die "Cannot open outfile4: $!\n";
print OUT4 "BarcodingTaxa";

foreach $samplePooledName (sort (keys (%pooledSamples))) {
	print OUT4 "\t".$samplePooledName;
}
print OUT4 "\n";

#FIRST PASS ONLY FOR GENUSBP >= 0.6
while ($parsed[$i]) {
	$line = $parsed[$i];
	chomp $line;
	@line = split(/\t/,$line);
	
	$sampleLine = $line[0];
	@sampleLine = split(/;/,$sampleLine);
	$site = $sampleLine[1];

	if (exists $samples{$site}) { ### only record genera from INV samples
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
		
#			if (exists $samples{$site}) {

			if ($genus eq $taxon || $genusNEW eq $taxon) {

				if ($genusBP >= $genusBPcutoff) {
			
					if ($size >= $minsizecutoff) {
			
						$oldvalue = $samples{$site}; #update samples counter for total Occurrence
						$newvalue = $oldvalue + 1;
						$samples{$site} = $newvalue;

						$samplePooledName = $sampleMappings{$site}; ### update pooled samples counter for total Occurrence
						$oldvalue2 = $pooledSamples{$samplePooledName};
						$newvalue2 = $oldvalue2 + 1; #increment
						$pooledSamples{$samplePooledName} = $newvalue2;

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

								$samplePooledName = $sampleMappings{$site}; ### update pooled samples counter for total Occurrence
								$oldvalue2 = $pooledSamples{$samplePooledName};
								$newvalue2 = $oldvalue2 + 1; #increment
								$pooledSamples{$samplePooledName} = $newvalue2;

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

	#totalOccurrencePerPooledSample, presenceAbsencePerPooledSample
	print OUT3 "$taxon";
	print OUT4 "$taxon";
	foreach $samplePooledName (sort (keys (%pooledSamples))) {
		$newvalue2 = $pooledSamples{$samplePooledName};
		print OUT3 "\t".$newvalue2;
		if ($newvalue2 > 0 ) {
			print OUT4 "\t1";
		}
		else {
			print OUT4 "\t0";
		}
		$pooledSamples{$samplePooledName} = "0"; #reset
	}
	print OUT3 "\n";
	print OUT4 "\n";

	$i=0;
}


close OUT;
close OUT2;
close OUT3;
close OUT4;

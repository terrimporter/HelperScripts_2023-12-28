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
my $seqSampleName;
my $fieldSampleName;
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
my %sampleMappings; #key = SeqSampleName, value = FieldSampleName
my %genera; #key = genus (with a BP >= 0.6), value = 1 to indicate presence
my %sampleTaxa; #key = sampleName, #key = taxon (genus>=0.6 or family_undef if genus<0.6 if genus not already present in sample), #value = 0 to increment
my %uniqueTaxa; #key = taxon, value = 1 to indicate presence

open (IN, "<", $ARGV[0]) || die "Cannot open infile2: $!\n";
@sampleMappings = <IN>;
close IN;

foreach $sampleLine (@sampleMappings) {
	chomp $sampleLine;
	@sampleLine = split(/\t/,$sampleLine);
	$seqSampleName = $sampleLine[0];
	$fieldSampleName = $sampleLine[1];
	$sampleMappings{$seqSampleName} = $fieldSampleName;
}

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@parsed = <IN2>;
close IN2;

open (OUT, ">>", "totalOccurrencePerSample.txt") || die "Cannot open outfile1: $!\n";
print OUT "BarcodingTaxa";

foreach $seqSampleName (sort(keys(%sampleMappings))) {
	chomp $seqSampleName;
	print OUT "\t".$seqSampleName;
}
print OUT "\n";

open (OUT2, ">>", "presenceAbsencePerSample.txt") || die "Cannot open outfile2: $!\n";
print OUT2 "BarcodingTaxa";

foreach $seqSampleName (sort(keys(%sampleMappings))) {
	chomp $seqSampleName;
	print OUT2 "\t".$seqSampleName;
}
print OUT2 "\n";

#BUILD HASHES

foreach $seqSampleName (sort(keys(%sampleMappings))) {
	chomp $seqSampleName;

	#FIRST round of parsing to build %genera for current sample
	while ($parsed[$i]) {
		$line = $parsed[$i];
		chomp $line;
		@line = split(/\t/,$line);

		$sampleLine = $line[0];
		@sampleLine = split(/;/,$sampleLine);
		$site = $sampleLine[1];

		#build %genera and %sampleTaxa for current sample
		if ($seqSampleName eq $site) {
			$genus = $line[14]; 
			$genusBP = $line[15]; 

			if ($genusBP >= $genusBPcutoff) {
				$genera{$genus} = "1"; #use to check against during SECOND round of parsing
				$sampleTaxa{$seqSampleName}{$genus} = "0"; #use to check against during THIRD round of parsing
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

	#SECOND round of parsing to build %sampleTaxa for all samples
	while ($parsed[$i]) {
		$line = $parsed[$i];
		chomp $line;
		@line = split(/\t/,$line);

		$sampleLine = $line[0];
		@sampleLine = split(/;/,$sampleLine);
		$site = $sampleLine[1];

		#continue to build %sampleTaxa for current sample
		if ($seqSampleName eq $site) {
			$genus = $line[14]; 
			$genusBP = $line[15]; 
			$family = $line[12];
			$familyBP = $line[13];

			if ($genusBP < $genusBPcutoff && $familyBP >= $familyBPcutoff) {

				if (exists $genera{$genus}) {
					$i++;
					next;
				}
				else {
					$taxon = $family."_undef";
					$sampleTaxa{$seqSampleName}{$taxon} = "0";
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
	%genera=(); #reset for next sample
}

#create hash of unique taxa only
foreach $seqSampleName (sort(keys(%sampleTaxa))) {
	foreach $taxon (keys %{$sampleTaxa{$seqSampleName}}) {
		if (exists $uniqueTaxa{$taxon}) {
			next;
		}
		else {
			$uniqueTaxa{$taxon} = 1;
		}
	}
}

#start incrementing taxon counters then print outfiles
foreach $taxon (sort(keys(%uniqueTaxa))) { ### considers all genera and family_undef per sample
	chomp $taxon;

	while ($parsed[$i]) {
		$line = $parsed[$i];
		chomp $line;
		@line = split(/\t/,$line);

		$sampleLine = $line[0];
		@sampleLine = split(/;/,$sampleLine);
		$site = $sampleLine[1];
		$size = $sampleLine[2];

		if (exists $sampleMappings{$site}) { ### only record genera & family_undef from INV samples

			if ($size >= $minsizecutoff) {

				$genus = $line[14];
				$genusBP = $line[15];

				if ($genusBP >= $genusBPcutoff) {

					if (exists $sampleTaxa{$site}{$genus}) {
						$oldvalue = $sampleTaxa{$site}{$genus};
						$newvalue = $oldvalue + 1;
						$sampleTaxa{$site}{$genus} = $newvalue;
					}
					else {
						print "site: $site genus:$genus does not exist1\n";
					}
				}
				else {
					$family = $line[12]; 	
					$familyBP = $line[13];
		
					if ($familyBP >= $familyBPcutoff) {
						$genusNEW = $family."_undef";
							
							if (exists $sampleTaxa{$site}{$genusNEW}) {
								$oldvalue = $sampleTaxa{$site}{$genusNEW};
								$newvalue = $oldvalue + 1;
								$sampleTaxa{$site}{$genusNEW} = $newvalue;
							}
							else {
								print "site: $site genus: $genusNEW does not exist probably because the genus already exists.\n";
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
	$i=0;
		
	#totalOccurrencePerSample
	print OUT "$taxon";
	#print values for each rep
	foreach $seqSampleName (sort(keys(%sampleTaxa))) {
		if (exists $sampleTaxa{$seqSampleName}{$taxon}) {
			$newvalue = $sampleTaxa{$seqSampleName}{$taxon};
			print OUT "\t".$newvalue;
		}
		else {
			print OUT "\t0";
			next;
		}
	}
	print OUT "\n";

	#presenceAbsencePerSample
	print OUT2 "$taxon";
	foreach $seqSampleName (sort(keys(%sampleTaxa))) {
		if (exists $sampleTaxa{$seqSampleName}{$taxon}) {
			if ($sampleTaxa{$seqSampleName}{$taxon} >= 1) {
				print OUT2 "\t1";
			}
			else {
				print OUT2 "\t0";
			}
		}
		else {
			print OUT2 "\t0";
		}
	}
	print OUT2 "\n";
}

close OUT;
close OUT2;

#!/usr/bin/perl
#Dec.8, 2016 by Terri Porter
#Script to make presence-absence and numOccurrence matrices from RDP classifier CO1v2 output
#Customize to process Blackbrook headers, genus assignments bp>0.6, else if family bp > 0.3 then family_undef
#USAGE perl make_matrices_from_RDP2.plx morphTaxonListYEAR.txt sampleNames.txt YEAR_rdp.out

#declare var
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
my $counter=0;
my $newcounter;
my $percentOccurrence;
#my $numReps;

#declare array
my @genusMorph;
my @sampleNames;
my @parsed;
my @line;
my @header;

#declare hash
%samples; #key = repSampleName ex. 1A, 1B, 1C, 3A, etc.; value="0" to be incremented/reset

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@genusMorph = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@sampleNames = <IN2>;
close IN2;

foreach $sampleName (@sampleNames) {
	chomp $sampleName;
	$samples{$sampleName} = "0"; 
}
#$numReps = scalar @sampleNames;

open (IN3, "<", $ARGV[2]) || die "Cannot open infile2: $!\n";
@parsed = <IN3>;
close IN3;

open (OUT, ">>", "totalOccurrencePerSample.txt") || die "Cannot open outfile1: $!\n";
print OUT "genusMorphList";

foreach $sampleName (@sampleNames) {
	chomp $sampleName;
	print OUT "\t".$sampleName;
}
print OUT "\n";

open (OUT2, ">>", "presenceAbsencePerSample.txt") || die "Cannot open outfile2: $!\n";
print OUT2 "genusMorphList";

foreach $sampleName (@sampleNames) {
	chomp $sampleName;
	print OUT2 "\t".$sampleName;
}
print OUT2 "\n";

#open (OUT3, ">>", "percentOccurrencePerTaxon.txt") || die "Cannot open outfile3: $!\n";
#print OUT3 "genusMorphList\tPercentOccurrence\n";

#open (OUT4, ">>", "presenceAbsencePerTaxon.txt" ) || die "Cannot open outfile4: $!\n";
#print OUT4 "genusMorphList\tPresenceAbsence\n";

foreach $genusMorph (@genusMorph) {
	chomp $genusMorph;
#	print "genusMorph: $genusMorph\n";#test

	while ($parsed[$i]) {
		$line = $parsed[$i];
#		print "line: $line\n";#test
		chomp $line;
		@line = split(/\t/,$line);

		$genus = $line[14]; ### CUSTOMIZE THIS HERE ### genus
#		print "genus: $genus\n"; #test
		$genusBP = $line[15]; ### CUSTOMIZE THIS HERE ### genusBP
#		print "bp: $bp\n"; #test

		if ($genusBP < $genusBPcutoff) {
			$family = $line[12]; 	
#			print "family: $family\n"; #test
			$familyBP = $line[13];
#			print "familyBP: $familyBP\n"; #test

			if ($familyBP >= $familyBPcutoff) {
				$genus = $family."_undef";
			}
		}

		$header = $line[0];
#		print "OTUheader: $header\n";#test
		@header = split(/;/,$header); ### CUSTOMIZE HERE ### Hearst & Blackbrook

		$site = $header[1];
#		print "site: $site\n";#test

		$size = $header[2];
#		print "size: $size\n";#test

		
		if ($genus eq $genusMorph) {

			#check for minimum bootstrap support cutoff value
			if ($genusBP >= $genusBPcutoff || $familyBP >= $familyBPcutoff) {
	
				#check for minimum OTU cluster size cutoff
				if ($size >= $minsizecutoff) {
			
					$oldvalue = $samples{$site};
					$newvalue = $oldvalue + 1; #increment
					$samples{$site} = $newvalue;

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
		}
		else {
			$i++;
			next;
		}

		$i++;
	}
		
	#totalOccurrencePerRep
	print OUT "$genusMorph";
	#print values for each rep
	foreach $sampleName (@sampleNames) {
		chomp $sampleName;
		print OUT "\t".$samples{$sampleName};
#		$reps{$repName} = "0"; #don't reset until last outfile processed
	}
	print OUT "\n";

	#presenceAbsencePerRep, percentOccurrencePerFamily, presenceAbsencePerFamily
	print OUT2 "$genusMorph";
#	print OUT3 "$genusMorph";
#	print OUT4 "$genusMorph";
	foreach $sampleName (@sampleNames) {
		chomp $sampleName;
		if ($samples{$sampleName} >= 1) {
			print OUT2 "\t1";
			$counter++;
		}
		else {
			print OUT2 "\t0";
		}
		$samples{$sampleName} = "0"; #reset
	}
	print OUT2 "\n";
#	$percentOccurrence = ($counter / $numReps) * 100;
#	print OUT3 "\t".$percentOccurrence."\n";
#	if ($counter > 0) {
#		print OUT4 "\t1\n";
#	}
#	else {
#		print OUT4 "\t0\n";
#	}
	$counter=0; #reset
#	$percentOccurrence=(); #clear

	$i=0;
}
close OUT;
close OUT2;
#close OUT3;
#close OTU4;

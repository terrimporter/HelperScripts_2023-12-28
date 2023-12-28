#!/usr/bin/perl
#Dec.8, 2016 by Terri Porter
#Script to make presence-absence and numOccurrence matrices from RDP classifier CO1v2 output
#Customize to process Hearst headers
#USAGE perl make_matrices_from_RDP.plx morphTaxonList.txt repNames.txt Hearst_CO1v2_concatenated.txt

#declare var
my $repName;
my $cutoff;
my $minsize;
my $taxonMorph;
my $line;
my $taxon; #family, genus, etc. 
my $OTUheader;
my $site; #field site replicate name
my $size; #OTUsize, i.e. number of reads clustered in OTU
my $bp; #RDP classifier bootstrap proportion (0-1.0)
my $oldvalue;
my $newvalue;
my $counter=0;
my $newcounter;
my $percentOccurrence;
my $numReps;

#declare array
my @taxonMorph;
my @repNames;
my @parsed;
my @line;

#declare hash
%reps; #key = repSampleName ex. 1A, 1B, 1C, 3A, etc.; value="0" to be incremented/reset

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@taxonMorph = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@repNames = <IN2>;
close IN2;

foreach $repName (@repNames) {
	chomp $repName;
	$reps{$repName} = "0"; 
}
$numReps = scalar @repNames;

open (IN3, "<", $ARGV[2]) || die "Cannot open infile2: $!\n";
@parsed = <IN3>;
close IN3;

print "Enter bootstrap support cutoff for genus rank (0-1.0, CO1v2 use 0.6 for 200bp fragment):\n";
$cutoff = <STDIN>;
chomp $cutoff;

print "Enter minimum OTU cluster size cutoff (ex. minsize=3 will exclude singletons and doubletons and retain a minimum cluster size of 3 reads):\n";
$minsize = <STDIN>;
chomp $minsize;

print "Will use a bootstrap support cutoff of $cutoff and a minimum OTU size of $minsize\n";

open (OUT, ">>", "totalOccurrencePerRep.txt") || die "Cannot open outfile1: $!\n";
print OUT "taxonMorphList";

foreach $repName (@repNames) {
	chomp $repName;
	print OUT "\t".$repName;
}
print OUT "\n";

open (OUT2, ">>", "presenceAbsencePerRep.txt") || die "Cannot open outfile2: $!\n";
print OUT2 "taxonMorphList";

foreach $repName (@repNames) {
	chomp $repName;
	print OUT2 "\t".$repName;
}
print OUT2 "\n";

open (OUT3, ">>", "percentOccurrencePerTaxon.txt") || die "Cannot open outfile3: $!\n";
print OUT3 "taxonMorphList\tPercentOccurrence\n";

open (OUT4, ">>", "presenceAbsencePerTaxon.txt" ) || die "Cannot open outfile4: $!\n";
print OUT4 "taxonMorphList\tPresenceAbsence\n";

foreach $taxonMorph (@taxonMorph) {
	chomp $taxonMorph;
#	print "taxonMorph: $taxonMorph\n";#test

	while ($parsed[$i]) {
		$line = $parsed[$i];
#		print "line: $line\n";#test
		chomp $line;
		@line = split(/\t/,$line);

		$taxon = $line[14]; ### CUSTOMIZE THIS HERE ### genus
#		print "taxon: $taxon\n"; #test

		$OTUheader = $line[0];
		print "OTUheader: $OTUheader\n";#test
		@OTUheader = split(/;/,$OTUheader); ### CUSTOMIZE HERE ### Hearst

		$site = $OTUheader[1];
		if ($site =~ /\d{1}\./) { ### customized for HEARST ###
			$site =~ /(\d{1})\./;
			$site = "DT0".$1;
		}
		print "site: $site\n";#test

		$size = $OTUheader[2];
		print "size: $size\n";#test

		$bp = $line[15]; ### CUSTOMIZE THIS HERE ### genusBP
		print "bp: $bp\n"; #test
	
		if ($taxon eq $taxonMorph) {

			#check for minimum bootstrap support cutoff value
			if ($bp >= $cutoff) {
	
				#check for minimum OTU cluster size cutoff
				if ($size >= $minsize) {
			
					$oldvalue = $reps{$site};
					$newvalue = $oldvalue + 1; #increment
					$reps{$site} = $newvalue;

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
	print OUT "$taxonMorph";
	#print values for each rep
	foreach $repName (@repNames) {
		chomp $repName;
		print OUT "\t".$reps{$repName};
#		$reps{$repName} = "0"; #don't reset until last outfile processed
	}
	print OUT "\n";

	#presenceAbsencePerRep, percentOccurrencePerFamily, presenceAbsencePerFamily
	print OUT2 "$taxonMorph";
	print OUT3 "$taxonMorph";
	print OUT4 "$taxonMorph";
	foreach $repName (@repNames) {
		chomp $repName;
		if ($reps{$repName} >= 1) {
			print OUT2 "\t1";
			$counter++;
		}
		else {
			print OUT2 "\t0";
		}
		$reps{$repName} = "0"; #reset
	}
	print OUT2 "\n";
	$percentOccurrence = ($counter / $numReps) * 100;
	print OUT3 "\t".$percentOccurrence."\n";
	if ($counter > 0) {
		print OUT4 "\t1\n";
	}
	else {
		print OUT4 "\t0\n";
	}
	$counter=0; #reset
	$percentOccurrence=(); #clear

	$i=0;
}
close OUT;
close OUT2;
close OUT3;

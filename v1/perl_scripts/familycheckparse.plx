#!/usr/bin/perl
#Dec.6, 2016 by Terri Porter
#Script to make presence-absence and numOccurrence matrices from parsed.txt from familycheck.plx
#USAGE perl familycheckparsed.plx familylist.txt repnames.txt parsed.txt

#declare var
my $repName;
my $cutoff;
my $minsize;
my $familyMorph;
my $line;
my $family; 
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
my @familyMorph;
my @repNames;
my @parsed;
my @line;

#declare hash
%reps; #key = repSampleName ex. 1A, 1B, 1C, 3A, etc.; value="0" to be incremented/reset

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@familyMorph = <IN>;
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

print "Enter bootstrap support cutoff for family rank (0-1.0, CO1v2 use 0.3 for 200bp fragment):\n";
$cutoff = <STDIN>;
chomp $cutoff;

print "Enter minimum OTU cluster size cutoff (ex. minsize=3 will exclude singletons and doubletons and retain a minimum cluster size of 3 reads):\n";
$minsize = <STDIN>;
chomp $minsize;

print "Will use a bootstrap support cutoff of $cutoff and a minimum OTU size of $minsize\n";

open (OUT, ">>", "totalOccurrencePerRep.txt") || die "Cannot open outfile1: $!\n";
print OUT "familyMorph";

foreach $repName (@repNames) {
	chomp $repName;
	print OUT "\t".$repName;
}
print OUT "\n";

open (OUT2, ">>", "presenceAbsencePerRep.txt") || die "Cannot open outfile2: $!\n";
print OUT2 "familyMorph";

foreach $repName (@repNames) {
	chomp $repName;
	print OUT2 "\t".$repName;
}
print OUT2 "\n";

open (OUT3, ">>", "percentOccurrencePerFamily.txt") || die "Cannot open outfile3: $!\n";
print OUT3 "familyMorph\tPercentOccurrence\n";

open (OUT4, ">>", "presenceAbsencePerFamily.txt" ) || die "Cannot open outfile4: $!\n";
print OUT4 "familyMorph\tPresenceAbsence\n";

foreach $familyMorph (@familyMorph) {
	chomp $familyMorph;

	while ($parsed[$i]) {
		$line = $parsed[$i];
		chomp $line;
		@line = split(/\t/,$line);

		$family = $line[0];
		$site = $line[1];
		$size = $line[2];
		$bp = $line[4];
	
		if ($family eq $familyMorph) {

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
	print OUT "$familyMorph";
	#print values for each rep
	foreach $repName (@repNames) {
		chomp $repName;
		print OUT "\t".$reps{$repName};
#		$reps{$repName} = "0"; #don't reset until last outfile processed
	}
	print OUT "\n";

	#presenceAbsencePerRep, percentOccurrencePerFamily, presenceAbsencePerFamily
	print OUT2 "$familyMorph";
	print OUT3 "$familyMorph";
	print OUT4 "$familyMorph";
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

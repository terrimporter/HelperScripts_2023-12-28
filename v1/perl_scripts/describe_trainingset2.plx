#!/usr/bin/perl
#need tabular summary of original lineages (not what they were misidentified as!)
#Nov. 17, 2016 no need to combine LOO partial results anymore, just describe the training set, for CO1 v2 need superkingdom rank added
#April 24, 2013 by Terri Porter
#Script to combine a bunch of RDP LOOCV output back into a single file
#usage perl combine_LOOCV.plx testNBC.taxonomy

use strict;
use warnings;

#declare var
my $dir;
my $i=0;
my $file;
my $filepath;
my $j=0;
my $line;
my $flag=0;
my $rank;
my $original;
my $newCorrect;
my $newTotal;
my $singleton=0;
my $rank2;
my $name;
my $tested;
my $misClassified;
my $originalRank;
my $originalTotal2;
my $originalTested;
my $originalMisClassified;
my $newTotal2;
my $newTested;
my $newMisClassified;
my $NAME;
my $correct95; #100-95
my $total95;
my $correct90; #94-90
my $total90;
my $correct80; #89-80
my $total80;
my $correct70; #79-70
my $total70;
my $correct60; #69-60
my $total60;
my $correct50; #59-50
my $total50;
my $correct40; #49-40
my $total40;
my $correct30; #39-30
my $total30;
my $correct20; #29-20
my $total20;
my $correct10; #19-10
my $total10;
my $correct0; #9-0
my $total0;
my $correct;
my $total;
my $originalCorrect95;
my $originalTotal95;
my $originalCorrect90;
my $originalTotal90;
my $originalCorrect80;
my $originalTotal80;
my $originalCorrect70;
my $originalTotal70;
my $originalCorrect60;
my $originalTotal60;
my $originalCorrect50;
my $originalTotal50;
my $originalCorrect40;
my $originalTotal40;
my $originalCorrect30;
my $originalTotal30;
my $originalCorrect20;
my $originalTotal20;
my $originalCorrect10;
my $originalTotal10;
my $originalCorrect0;
my $originalTotal0;
my $originalCorrect;
my $originalTotal;
my $newCorrect95;
my $newTotal95;
my $newCorrect90;
my $newTotal90;
my $newCorrect80;
my $newTotal80;
my $newCorrect70;
my $newTotal70;
my $newCorrect60;
my $newTotal60;
my $newCorrect50;
my $newTotal50;
my $newCorrect40;
my $newTotal40;
my $newCorrect30;
my $newTotal30;
my $newCorrect20;
my $newTotal20;
my $newCorrect10;
my $newTotal10;
my $newCorrect0;
my $newTotal0;
my $numPartFiles=1; ### edit this ### OMIT THIS
my $adjTotal;
my $misclassified=0;
my $misclassifiedAcc;
my $superkingdom; #new
my $kingdom; #for each misclassified seq, keep the erroneous classification
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $superkingdomVal;
my $kingdomVal;
my $phylumVal;
my $classVal;
my $orderVal;
my $familyVal;
my $genusVal;
my $string;
my $new;
my $superkingdomO;
my $kingdomO; #for each misclassified seq, retain what the correct classification should have been
my $phylumO;
my $classO;
my $orderO;
my $familyO;
my $genusO;
my $stringO;

#declare array
my @files;
my @in;
my @line;
my @original;
my @tax;
my @rank = ('superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus');
my @string;
my @stringO;

#declare hash
my %one; #indexed by rank
my %table; #indexed by name
my %misclassifiedOriginal; #for each misclassified seq, record original lineage for table 2 recreation 90+
my %misclassified; #indexed by accession, for each misclassified seq, record misclassified lineage for misclassified section
my %table2; #greater than 0.90 including singletons
my %table3; #any bs including singletons

print "Print name of directory containing LOO file testNBC.out including final /:\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Cannot open directory: $!\n";
@files = readdir(DIR);
closedir(DIR);

open (TESTMIS, ">>", "loo_misclassified_singletons.txt") || die "Error cannot open loo_misclassified.txt:$!\n";

#open (TESTSINGLE, ">>", "test_singletons.txt") || die "Error cannot ope test_singletons.txt: $!\n";

while ($files[$i]) {
	$file = $files[$i];
	chomp $file;
#	print "processing $file\n";

	if ($file !~ /^\./) {
		$filepath = $dir.$file;
		print "Processing $filepath\n";
	}
	else {
		$i++;
		next;
	}

	open (IN, "<", $filepath) || die "Cannot open file: $!\n";
	@in = <IN>;
	close IN;

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($flag == 0) {
			if ($line =~ /1: number of correct assigned sequences/) {
				$flag = 1;
				print "Found first table\n";
				$j++;
				next;
			}
			else {
				$j++;
				next;
			}
		}

		if ($flag == 1) {
#			print "got to flag 1\n";
			if ($line =~ /^Level/) {
#				print "$line\n";
				$j++;
				next;
			}
			elsif ($line =~ /^(superkingdom|phylum|kingdom|order|family|class|genus)/) {
#				print "Processing first table\n";
				@line = split(' ', $line);
				$rank = $line[0]; #rank
				$correct95 = $line[1];
				$total95 = $line[2];
				$correct90 = $line[3];
				$total90 = $line[4];
				$correct80 = $line[5];
				$total80 = $line[6];
				$correct70 = $line[7];
				$total70 = $line[8];
				$correct60 = $line[9];
				$total60 = $line[10];
				$correct50 = $line[11];
				$total50 = $line[12];
				$correct40 = $line[13];
				$total40 = $line[14];
				$correct30 = $line[15];
				$total30 = $line[16];
				$correct20 = $line[17];
				$total20 = $line[18];
				$correct10 = $line[19];
				$total10 = $line[20];
				$correct0 = $line[21];
				$total0 = $line[22];
				$correct = $line[23];
				$total = $line[24];

				if (exists $one{$rank}) {
#					print "FOUND EXISTING ONE RANKS\n";
					$original = $one{$rank};
					@original = split(/\|/, $original);

					$originalCorrect95 = $original[0];
					$originalTotal95 = $original[1];
					$originalCorrect90 = $original[2];
					$originalTotal90 = $original[3];
					$originalCorrect80 = $original[4];
					$originalTotal80 = $original[5];
					$originalCorrect70 = $original[6];
					$originalTotal70 = $original[7];
					$originalCorrect60 = $original[8];
					$originalTotal60 = $original[9];
					$originalCorrect50 = $original[10];
					$originalTotal50 = $original[11];
					$originalCorrect40 = $original[12];
					$originalTotal40 = $original[13];
					$originalCorrect30 = $original[14];
					$originalTotal30 = $original[15];
					$originalCorrect20 = $original[16];
					$originalTotal20 = $original[17];
					$originalCorrect10 = $original[18];
					$originalTotal10 = $original[19];
					$originalCorrect0 = $original[20];
					$originalTotal0 = $original[21];
					$originalCorrect = $original[22];
					$originalTotal = $original[23];

					$newCorrect95 = $originalCorrect95+$correct95;
					$newTotal95 = $originalTotal95+$total95;
					$newCorrect90 = $originalCorrect90+$correct90;
					$newTotal90 = $originalTotal90+$total90;
					$newCorrect80 = $originalCorrect80+$correct80;
					$newTotal80 = $originalTotal80+$total80;
					$newCorrect70 = $originalCorrect70+$correct70;
					$newTotal70 = $originalTotal70+$total70;
					$newCorrect60 = $originalCorrect60 + $correct60;
					$newTotal60 = $originalTotal60 + $total60;
					$newCorrect50 = $originalCorrect50 + $correct50;
					$newTotal50 = $originalTotal50 + $total50;
					$newCorrect40 = $originalCorrect40 + $correct40;
					$newTotal40 = $originalTotal40 + $total40;
					$newCorrect30 = $originalCorrect30 + $correct30;
					$newTotal30 = $originalTotal30 + $total30;
					$newCorrect20 = $originalCorrect20 + $correct20;
					$newTotal20 = $originalTotal20 + $total20;
					$newCorrect10 = $originalCorrect10 + $correct10;
					$newTotal10 = $originalTotal10 + $total10;
					$newCorrect0 = $originalCorrect0 + $correct0;
					$newTotal0 = $originalTotal0 + $total0;
					$newCorrect = $originalCorrect+$correct;
					$newTotal = $originalTotal+$total;
					$one{$rank} = 	$newCorrect95."|".$newTotal95."|".
								  	$newCorrect90."|".$newTotal90."|".
									$newCorrect80."|".$newTotal80."|".
									$newCorrect70."|".$newTotal70."|".
									$newCorrect60."|".$newTotal60."|".
									$newCorrect50."|".$newTotal50."|".
									$newCorrect40."|".$newTotal40."|".
									$newCorrect30."|".$newTotal30."|".
									$newCorrect20."|".$newTotal20."|".
									$newCorrect10."|".$newTotal10."|".
									$newCorrect0."|".$newTotal0."|".
									$newCorrect."|".$newTotal;
				}
				else {
					$one{$rank} = 	$correct95."|".$total95."|".
									$correct90."|".$total90."|".
									$correct80."|".$total80."|".
									$correct70."|".$total70."|".
									$correct60."|".$total60."|".
									$correct50."|".$total50."|".
									$correct40."|".$total40."|".
									$correct30."|".$total30."|".
									$correct20."|".$total20."|".
									$correct10."|".$total10."|".
									$correct0."|".$total0."|".
									$correct."|".$total;
				}
				$j++;
				next;
			}
			else {
#				print "$line\n";
				$flag=6;
				$j++;
				next;
			}
		}

		if ($flag == 6) {
#			print "got to flag 6\n";
			if ($line =~ /\*\*misclassified sequences:/) { #account for typo in NBC output, corrected in v2.6
#				print "found misclassified sequence section\n"; #test
				$flag=7;
				$j++;
				next;
			}
			else {
				$j++;
				next;
			}
		}
		
		if ($flag == 7 || $flag == 3) {
#			print "got to flag 7\n";
			if ($line =~ /^SEQ:/) {
				print TESTMIS $line."\n";
#				print "Processing misclassified sequences\n";
				$line =~ /SEQ:\s{1}(\w+)/;
				$misclassifiedAcc = $1;
				$j++;

				$line = $in[$j]; #correct annotation
				chomp $line;
				print TESTMIS $line."\n";
				@line = split(' ',$line); 
				$superkingdomO = $line[0];
				$kingdomO = $line[1];
				$phylumO = $line[2];
				$classO = $line[3];
				$orderO = $line[4];
				$familyO = $line[5];
				$genusO = $line[6];
				$j++;
				
				$line = $in[$j]; #incorrect annotation
				chomp $line;
				print TESTMIS $line."\n";
				@line = split(' ',$line); 
				$superkingdom = $line[0];
				$kingdom = $line[1];
				$phylum = $line[2];
				$class = $line[3];
				$order = $line[4];
				$family = $line[5];
				$genus = $line[6];
				$j++;
				
				$line = $in[$j]; #assignment support
				chomp $line;
				print TESTMIS $line."\n";
				@line = split(' ',$line);
				$superkingdomVal = $line[0];
				$kingdomVal = $line[1];
				$phylumVal = $line[2];
				$classVal = $line[3];
				$orderVal = $line[4];
				$familyVal = $line[5];
				$genusVal = $line[6];
	
				$stringO = $misclassifiedAcc."\t\t".$superkingdomO."\tsuperkingdom\t".$superkingdomVal."\t".$kingdomO."\tkingdom\t".$kingdomVal."\t".$phylumO."\tphylum\t".$phylumVal."\t".$classO."\tclass\t".$classVal."\t".$orderO."\torder\t".$orderVal."\t".$familyO."\tfamily\t".$familyVal."\t".$genusO."\tgenus\t".$genusVal;
				if (exists $misclassifiedOriginal{$misclassifiedAcc}) {
#					print "found $misclassifiedAcc more than once in detail sections for different ranks\n"; #only process once!
				}
				else {
					$misclassifiedOriginal{$misclassifiedAcc} = $stringO;
				}

				$string = $misclassifiedAcc."\t\t".$superkingdom."\tsuperkingdom\t".$superkingdomVal."\t".$kingdom."\tkingdom\t".$kingdomVal."\t".$phylum."\tphylum\t".$phylumVal."\t".$class."\tclass\t".$classVal."\t".$order."\torder\t".$orderVal."\t".$family."\tfamily\t".$familyVal."\t".$genus."\tgenus\t".$genusVal;

				if (exists $misclassified{$misclassifiedAcc}) {
#					print "found $misclassifiedAcc more than once in detail sections for different ranks\n"; #only process once!
				}
				else {
					$misclassified{$misclassifiedAcc} = $string;
				}

				$j++;
				next;
			}

			else {
#				print TESTMIS $line."\n";
				$flag = 2;
				$j++;
				next;
			}
		}

		if ($flag == 2) {
#			print "got to flag 2\n";
			if ($line =~ /\*\*singleton sequences:/) {
				print "Found singleton sequences section\n";
				$flag=3;
				$j++;
				next;
			}
			elsif ($line =~ /\*\*misclassified sequences group by taxon/) { ### new in v2.6, was ^Rank v2.5
#				print "Found second table.\n"; #test
				$flag=5;
				$j+=2;
				next;
			}
			else {
				$j++;
				next;
			}
		}

		if ($flag == 5) {
			if ($line =~ /^\S+/) {
				@line = split(' ',$line);
				$rank2 = $line[0];
				$name = $line[1];
				$total = $line[2];
				$tested = $line[3];
				$misClassified = $line[4];
#				print "name:$name\n"; #test

				if (exists $table{$name}) {
#					print "FOUND EXISTING NAMES\n"; #test
					$original = $table{$name};
					@original = split(/\|/,$original);
					$originalRank = $original[0];
					$originalTotal2 = $original[1];
					$originalTested = $original[2];
					$originalMisClassified = $original[3];
					$newTotal2 = $originalTotal2+$total;
					$newTested = $originalTested+$tested;
					$newMisClassified = $originalMisClassified+$misClassified;
					$table{$name} = $rank2."|".$newTotal2."|".$newTested."|".$newMisClassified;
				}
#				elsif ($line =~ /^\s+/) {
#					$flag=4; ## do not parse new ROC or Area under the Curve sections new in v2.6
#					$j++;
#				}
				else {
					$table{$name} = $rank2."|".$total."|".$tested."|".$misClassified;
				}
				@line = ();
				$rank2=();
				$name=();
				$total=();
				$tested=();
				$misClassified=();
				$original=();
				@original=();
				$originalRank=();
				$originalTotal2=();
				$originalTested=();
				$originalMisClassified=();
				$newTotal2=();
				$newTested=();
				$newMisClassified=();
				$j++;
				next;
			}
			else {
				$flag = 4; ## do not parse any further
				print "Found end of second table\n";
				$j++;
			}
		}

		if ($flag == 4) {
			$j++; #flip through end of file but do nothing, new for v2.6
		}
	}
	$j=0;
	$i++;
	$flag=0;
	@in=();
}
$i=0;
$j=0;

close TESTMIS;
#close TESTSINGLE;

open (OUT, ">>", "loo.txt") || die "Error cannot open outfile: $!\n";

print OUT "There were $misclassified misclassified genera printed to outfile loo.table.\n\n";

open (OUT2,">>","loo_misclassified.table") || die "Error cannot open outfile2:$!\n";

while ( ($misclassifiedAcc,$string) = each(%misclassifiedOriginal) ) { ### NEW show original classificiation, not what it was misclassified as!!! ##
	print OUT2 $string."\n";
}
close OUT2;

print OUT "There are $singleton singleton genera not used in the following analyses.\n\n";

print OUT "Statistics for each hierarchy level:\n";
print OUT "Rank\tCorrect 100-95\tTotal 100-95\tCorrect 94-90\tTotal 94-90\tCorrect 89-80\tTotal 89-80\tCorrect 79-70\tTotal 79-70\tCorrect 69-60\tTotal 69-60\tCorrect 59-50\tTotal 59-50\tCorrect 49-40\tTotal 49-40\tCorrect 39-30\tTotal 39-30\tCorrect 29-20\tTotal 29-20\tCorrect 19-10\tTotal 19-10\tCorrect 9-0\tTotal 9-0\tCorrect\tTotal\n";

while ($rank[$i]) {
	$rank = $rank[$i];
	chomp $rank;
#	while( ($rank,$line) = each (%one) ) {
	if (exists $one{$rank}) {
		$line = $one{$rank};
		print OUT "$rank\t";
		$line =~ s/\|/\t/g;
		print OUT "$line\n";
#		@line = split(/\|/,$line);
#		$correct95 = $line[0];
#		$total95 = $line[1];
#		$correct90 = $line[2];
#		$total90 = $line[3];
#		$correct = $line[4];
#		$total = $line[5];
#		print OUT "$correct95\t$total95\t$correct90\t$total90\t$correct\t$total\n";
	}
	else {
		print "Error cannot find data for $rank in loo in Table 1\n";
	}
	$i++;
}
$i=0;

print OUT "\nStatistics for each taxon in training set:\n";
print OUT "Rank\tName\tTotal\tTested\tMisClassified\n";

#print this out in the same order as the testNBC.taxonomy file

open (TAX, "<", $ARGV[0]) || die "Error cannot open taxonomy file: $!\n";
@tax = <TAX>;
close TAX;

while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	@line = split(/\*/, $line);
	$NAME = $line[1];
#	print "NAME:$NAME\n";

#while ( ($name,$line) = each (%table) ) {
	if (exists $table{$NAME} ) {
		$line = $table{$NAME};
		@line = split(/\|/, $line);
		$rank = $line[0];
		$total = $line[1];
		$adjTotal = $total/$numPartFiles;
		$tested = $line[2];
		$misClassified = $line[3];
		print OUT "$rank\t$NAME\t$adjTotal\t$tested\t$misClassified\n";
	}
	else {
		print "Couldn't find data for $NAME in any loo files\n";
	}
	$i++;
}
$i=0;

#recreate table 2 for 90+ only at order rank only
while ( ($misclassifiedAcc,$stringO) = each(%misclassifiedOriginal) ) {
	@stringO = split(' ',$stringO);
	$kingdomO = $stringO[1];
	$phylumO = $stringO[4];
	$classO = $stringO[7];
	$orderO = $stringO[10];
	$familyO = $stringO[13];
	$genusO = $stringO[16];

	$kingdomVal = $stringO[3];
	$phylumVal = $stringO[6];
	$classVal = $stringO[9];
	$orderVal = $stringO[12];
	$familyVal = $stringO[15];
	$genusVal = $stringO[18];

	#use this as table 2 instead of default output, because this will include singletons
	if (exists $table3{$orderO}) {
		$original = $table3{$orderO};
		$new = $original+1;
		$table3{$orderO} = $new;
	}
	else {
		$table3{$orderO} = 1;
	}

	#use this for filtered 90% bs output for table two, that also includes singletons
	if ($genusVal >= 0.90) {
		if (exists $table2{$orderO}) {
			$original = $table2{$orderO};
			$new = $original+1;
			$table2{$orderO} = $new;
		}
		else  {
			$table2{$orderO} = 1;
		}
	}
}

open (OUT3, ">>", "loo_misclassified_order_table2_90_countsingletons.txt") || die "Error cannot open third outfile: $!\n";

while ( ($orderO,$stringO) = each(%table2) ) {
	print OUT3 "$orderO\t$stringO\n";
}
close OUT3;

open (OUT4, ">>", "loo_misclassified_order_table2_countsingletons.txt") || die "Error cannot open fourth outfile: $!\n";

while (($orderO,$stringO) = each(%table3)) {
	print OUT4 "$orderO\t$stringO\n";
}

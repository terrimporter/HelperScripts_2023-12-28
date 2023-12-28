#!/usr/bin/perl
#Jun. 30/17 by Terri Porter
#Script to calculate false positives from merged blast out files with query coverage already calculated
#USAGE perl find_false_positives_BLAST.plx file_merged_qcov.txt mytrainseq.fasta

use strict;
use warnings;

#declare var
my $primer;
my $i=0;
my $line;
my $totalNumQueries;
my $queryID;
my $pident;
my $qcov;
my $QueryGenus;
my $RefID;
my $RefGenus;
my $falsepositive=0; #incorrect with high blast stats
my $truepositive=0; #correct with high blast stats
my $falsenegative=0; #correct with low blast stats
my $truenegative=0; #incorrect with low blast stats

#declare array
my @in;
my @filenamesplit;
my @split;
my @ref;

#declare hash
my %ref; #key=RefID, value=RefGenus

open (IN, "$ARGV[0]") || die "Error cannot open infile: $!\n";
@in=<IN>;
close IN;

@filenamesplit=split(/_/,$ARGV[0]);
$primer=shift @filenamesplit;

$totalNumQueries=scalar(@in);

open (REF,"$ARGV[1]") || die "Error cannot open ref infile: $!\n";
@ref=<REF>;
close REF;

#hash reference taxonomy
while ($ref[$i]) {
	$line=$ref[$i];
	chomp $line;

	if ($line=~/^>/) {
		@split=split(/;/,$line);
		$RefID=shift @split;
		$RefID=~s/>//g;
		$RefID=~s/\s+\w+$//g;
#		print "RefID:$RefID\t"; #test
		$RefGenus=pop @split;
#		print "RefGenus:$RefGenus\n"; #test
		$ref{$RefID}=$RefGenus;
	}
	else {
		$i++;
		next;
	}
	$i++;
	@split=();
}
$i=0;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@split=split(/\|/,$line);
	$queryID=$split[0];
	$pident=$split[7];
	$qcov=$split[27];
	$QueryGenus=$split[24];

#	print $queryID."\t".$pident."\t".$qcov."\n"; #test

	if ($pident >=95 && $qcov >= 85) {
		if (exists $ref{$queryID}){
			$RefGenus=$ref{$queryID};
#			print "RefGenus:$RefGenus\n"; #test
#			print "found queryid in hash\n"; #test
			
			if ($QueryGenus eq $RefGenus) {
				$truepositive+=1;
				$i++;
				next;
			}
			else {
				$falsepositive+=1;
			}
		}
		else {
			print "Error cannot find queryID $queryID in ref hash\n";
		}
#		print $QueryGenus."\n"; #test
	}
	else {
		if (exists $ref{$queryID}) {
			$RefGenus=$ref{$queryID};
			if ($QueryGenus eq $RefGenus) {
				$falsenegative+=1;
				$i++;
				next;
			}
			else {
				$truenegative+=1;
			}
		}
		else {
			print "Error cannot find queryID $queryID in ref hash\n";
		}
	}
$i++;
}
$i=0;

print "$primer\t$totalNumQueries\t$truepositive\t$falsenegative\t$truenegative\t$falsepositive\n";

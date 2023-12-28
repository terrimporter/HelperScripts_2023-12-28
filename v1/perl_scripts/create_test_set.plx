#!/usr/bin/perl
#June 6, 2012 by Terri Porter
#Parse through a set of alignment files, and just keep the ones with >= 400 aligned baits
#usage perl create_test_set.plx

use strict;
use warnings;

#declare var
my $file;
my $i=0;
my $j=0;
my $line;
my $output;
my $k;
my $gb;
my $line2;
my $alnBait;
my $ref;
my $minBaits = 400; ### Edit minimum number of aligned baits required to be in file before parsing
my $range;
my $sampleSize = 1; ### Edit sample size here
my $outfile;
my $randIndex;
my $random;

#declare array
my @output;
my @in;
my @file;
my @keys;
my @lines;

#declare hash
my %aln;
my %lines;

@output = qx(ls | grep 'fasta\$');

#parse through each file individually
while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	@file = split(/\./, $file);
	$ref = $file[0];

	$output = qx(grep '>' $file | wc -l);

	if ($output >= $minBaits) {
		
		open (IN,"<",$file) || die "Error cannot open $file: $!\n";
		@in = <IN>;
		close IN;
		
		#add baits to hash
		while ($in[$j]) {
			$line = $in[$j];
			chomp $line;

			if ($line =~ /^>/) {
				$line =~ s/^>//;
				$gb = $line;

				$k = $j+1;
				$line2 = $in[$k];
				chomp $line2;

				if ($line2 =~ /\S/){
					$alnBait = $line2;
					$aln{$gb} = $alnBait;
					$j+=2;
				}
			}
			else {
				$j++;
			}
			$line=();
			$gb=();
			$k=();
			$line2=();
			$alnBait=();
		}
		$j=0;
		
		#add gb keys to array for ramdom sampling below
		@keys = keys(%aln);
		$range = scalar(@keys);

		#randomly sample without replacement scheme, use hash to ensure unique random keys
		while (scalar(keys %lines) < $sampleSize) {
			$random = int(rand($range)); #zero (included) to top end (not inclusive)
			#$random = $random+1; #one (included) to top end (included)
			$lines{$random} = 1; #ensures that each key is unique!
		}

		#add unique random keys to array to guide random gb key sampling
		@lines = keys(%lines);

		$outfile = $ref.".sim";
		open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

		#proceed to randomly sample gb keys
		while ($lines[$j]) {
			$randIndex = $lines[$j];
			$gb = $keys[$randIndex];#randomly sample from array of gb keys

			if (exists $aln{$gb}) {
				$alnBait = $aln{$gb};
				print OUT ">$gb\n$alnBait\n";
			}
			else {
				print "Error in file $file, problem with $gb\n";
			}

			$j++;
			$randIndex=();
			$gb=();
			$alnBait=();
			$outfile=();
		}
		$j=0;
		close OUT;
	}

	$i++;
	$file=();
	@file=();
	$ref=();
	$output=();
	@in=();
	@keys=();
	$range=();
	%lines=();
	$random=();
	@lines=();
	%aln=();
}
$i=0;



#!/usr/bin/perl
#Sept. 16, 2011 edit to work with features.txt from ITS_entrez_grab_gb.plx
#Sept. 15, 2011 by Terri Porter
#Script to parse features.txt from LSU_entrez_grab_gb.plx to prep file for [R] frequency histogram of sequence lengths
#usage perl parse_features.plx features.txt

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $i=0;
my $line;
my $organism;
my $spacer1;
my $ITSseq;
my $spacer2;
my $ITS;
my $length_spacer1;
my $length_ITSseq;
my $length_spacer2;
my $length;
my $num;
my $min;
my $max;
my $mean;

#declare array
my @features;
my @line;
my @spacer1;
my @ITSseq;
my @spacer2;
my @ITS;
my @length;

open (IN,"<",$ARGV[0]) || die "Error cannot open features.txt: $!\n";
@features = <IN>;
close IN;

open (OUT,">>","length.txt") || die "Error cannot write to length.txt: $!\n";

while ($features[$i]) {
	$line = $features[$i];
	chomp $line;
	
	if ($line =~ /^\S+/) {
		@line = split(/\t/,$line);
		$organism = $line[1];
		if ($organism !~ /(ID\tOrganism\tITS sequence|sp\.|cf\.|aff\.|uncultured|endophyte)/) {
			$spacer1 = $line[2];
			$ITSseq = $line[3];
			$spacer2 = $line[4];

			if ($spacer1 !~ /nil/) {
				if ($spacer1 =~ /\S+/) {
					@spacer1 = split(//,$spacer1);
					$length_spacer1 = scalar(@spacer1);
				}
			}
			elsif ($spacer1 =~ /nil/) {
				$length_spacer1 = 0;
				$spacer1=();
			}
			if ($ITSseq !~ /nil/) {
				if ($ITSseq =~ /\S+/) {
					@ITSseq = split(//,$ITSseq);
					$length_ITSseq = scalar(@ITSseq);
				}
			}
			elsif ($ITSseq =~ /nil/) {
				$length_ITSseq = 0;
				$ITSseq=();
			}
			if ($spacer2 !~ /nil/) {
				if ($spacer2 =~ /\S+/) {
					@spacer2 = split(//,$spacer2);
					$length_spacer2 = scalar(@spacer2);
				}
			}
			elsif ($spacer2 =~ /nil/) {
				$length_spacer2 = 0;
				$spacer2=();
			}
			if ($spacer1) {
				if ($ITSseq) {
					if ($spacer2) {
						$ITS = $spacer1.$ITSseq.$spacer2;
					}
					else {
						$ITS = $spacer1.$ITSseq;
					}
				}
				else {
					if ($spacer2) {
						$ITS = $spacer1.$spacer2;
					}
					else {
						$ITS = $spacer1;
					}
				}
			}
			else {
				if ($ITSseq) {
					if ($spacer2) {
						$ITS=$ITSseq.$spacer2;
					}
					else {
						$ITS = $ITSseq;
					}
				}
				else {
					if ($spacer2) {
						$ITS=$spacer2;
					}
					else {
						$ITS="nil";
					}
				}
			}
			if ($ITS !~/nil/) {
				@ITS = split(//,$ITS);
				$length = scalar(@ITS);
				if ($length > 0) {
					print OUT "$length\n";
					push(@length,$length);
				}
				elsif ($length ==0) {
					$i++;

					$spacer1=();
					$ITSseq=();
					$spacer2=();
					@spacer1=();
					@ITSseq=();
					@spacer2=();
					@ITS=();
					@line=();
					$length_spacer1=();
					$length_ITSseq=();
					$length_spacer2=();
					$length=();
					next;
				}
			}
		}
	}
	$i++;
	
	$spacer1=();
	$ITSseq=();
	$spacer2=();
	@spacer1=();
	@ITSseq=();
	@spacer2=();
	@ITS=();
	@line=();
	$length_spacer1=();
	$length_ITSseq=();
	$length_spacer2=();
	$length=();
}

$num = scalar(@length);
$min = min(@length);
$max = max(@length);
$mean = mean(@length);

print "Num seqs = $num\nMin length = $min\nMax length = $max\nMean length = $mean\n";

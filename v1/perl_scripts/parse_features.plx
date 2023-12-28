#!/usr/bin/perl
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
my $seq;
my $length;
my $num;
my $min;
my $max;
my $mean;

#declare array
my @features;
my @line;
my @seq;
my @length;

open (IN,"<",$ARGV[0]) || die "Error cannot open features.txt: $!\n";
@features = <IN>;
close IN;

open (OUT,">>","length.txt") || die "Error cannot write to length.txt: $!\n";

while ($features[$i]) {
	$line = $features[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$organism = $line[1];

	if ($organism !~ /(ID\tOrganism\tLSU sequence|sp\.|cf\.|aff\.|uncultured|endophyte)/) {
	
		$seq = $line[2];

		if ($seq !~ /nil/) {
			if ($seq =~ /\S+/) {
				@seq = split(//,$seq);
				$length = scalar(@seq);
				print OUT "$length\n";
				push(@length,$length);
			}
		}
	}
	$i++;

	@line=();
	$seq=();
	@seq=();
	$length=();
}

$num = scalar(@length);
$min = min(@length);
$max = max(@length);
$mean = mean(@length);

print "Num seqs = $num\nMin length = $min\nMax length = $max\nMean length = $mean\n";

#!/usr/bin/perl
#Sept. 12, 2013 edited to retain size of each cluster
#Jan.9, 2013 by Terri Porter
#Script to strip rank and xy coord from readid names from .blastn files (originally found in .fas files from Joel & Nagissa)
#usage perl parse_blastn_queryid.plx file.blastn

use strict;
use warnings;

#declare var
my $filename;
my $i=0;
my $line;
my $queryline;
my $queryid_old;
my $queryid_new;
my $size; #NEW

#declare array
my @in;
my @queryline;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0].".parsed";
open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";


while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /Query=/) {
		$queryline = $line;
		$queryline =~ s/\s+//g; #remove spaces
		$queryline =~ s/-/_/g; #replace dashes with underscores
		$queryline =~ s/=/_/g; #replace = with underscore
#		@queryline = split(/=/,$queryline); #split on the "="
#		$queryid_old = $queryline[1];
#		$size = $queryline[2];
#		$size =~ s/;//g;
#		if ($queryid_old =~ /^\S{14}/) {
#			$queryid_old =~ /^(\S{14})/;
#			$queryid_new = $1;
			print OUT "Query= $queryline\n";
			#print "$queryid_old\t$queryid_new\n";
#			$queryid_new=();
#		}
	}
	else {
		print OUT "$line\n";
	}
	$queryline=();
	@queryline=();
	$queryid_old=();
	$size=();

	$i++;
	$line=();
}
$i=0;

close OUT;

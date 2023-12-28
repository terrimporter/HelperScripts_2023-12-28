#!/usr/bin/perl
#Feb. 15, 2013 by Terri Porter
#Script to reformat taxon labels from raxml so they work for phylip
#usage perl reformat_raxml_bootstrap_for_phylip.plx RAxML_bootstrap.OUT

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $j=0;
my $part;
my $genus;
my $species;
my $gb;
my $count_start;
my $count_end;
my $label;
my $new;
my $newline;

#declare array
my @in;
my @line;
my @part;
my @newline;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "RAxML_bootstrap.NEW") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/,/,$line);

	while ($line[$j]) {
		$part = $line[$j];
		@part = split(/\_/,$part);
		
		
		if ($part[0]) {
			$genus = $part[0];
			if ($genus =~ /^\(/) {
				$count_start = () = $genus =~ /\(/g; #perl goatse operator
#				print "opening brackets: $count_start\n"; #test
				$genus =~ s/\(//g; #remove starting brackets
			}
			else {
				$count_start=0;
			}
		}
		if ($part[1]) {
			$species = $part[1];
		}
		if ($part[2]) {
			$gb = $part[2];
			if ($gb =~ /\)$/) {
				$count_end = () = $gb =~ /\)/g; #perl goatse operator
#				print "closing brackets: $count_end\n"; #test
				$gb =~ s/\)//g; #remove ending brackets
			}
			else {
				$count_end=0;
			}
		}

		$label = $gb."_".$genus."_".$species;
		$new = ("(" x $count_start).$label.(")" x $count_end);
#		print "newline $new\n";
		push(@newline, $new);

		$part=();
		@part=();
		$genus=();
		$gb=();
		$count_start=();
		$count_end=();
		$new=();

		$j++;
	}
	$j=0;

	$newline = join(',',@newline);
	print OUT "$newline;\n";
	@newline=();
	$newline=();

	$i++;

	$line=();
	@line=();
}
$i=0;

close OUT;	

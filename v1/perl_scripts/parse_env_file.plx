#!/usr/bin/perl
#Nov. 19, 2013 by Terri Porter
#Script to change blocks to subsites ex. PAD1A to PAD1 in env.txt file
#usage perl parse_env_file.plx env.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $otuid;
my $sample;
my $abund;

#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "env.txt.noblock") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /size/) {
		@line = split(/\t/,$line);
		$otuid = $line[0];
		$sample = $line[1];
		$abund = $line[2];

		$sample =~ s/(A|B|C)$//g; #remove block designation

		print OUT "$otuid\t$sample\t$abund\n";
	}
	$i++;
	$line=();
	@line=();
	$otuid=();
	$sample=();
	$abund=();
}
$i=0;
close OUT;

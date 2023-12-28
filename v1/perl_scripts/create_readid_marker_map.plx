#!/usr/bin/perl
#Oct.7, 2013 by Terri Porter
#Script to create readid_marker.map (for heatmaps)
#usage perl create_readid_marker_map.plx readid_sample.map

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $readid;
my $marker="rbcLF"; ### enter this manually, and customize for each marker ###

#declare array
my @in;
my @line;

open (IN, "<", $ARGV[0]) || die "Error cannot open readid_sample.map: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "readid_marker.map") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if (length($line)>0){
		@line = split(/\t/,$line);
		$readid = $line[0];
		print OUT "$readid\t$marker\n";
	}
	$i++;
	$line=();
	@line=();
	$readid=();
}
$i=0;
close OUT;

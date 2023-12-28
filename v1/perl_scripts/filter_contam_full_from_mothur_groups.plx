#!/usr/bin/perl
#June 12, 2012 by Terri Porter
#Script to filter contam.full (from filter_contam_from_mothur_list.plx) from file.groups (from make_mothur_groups_file.plx)
#usage perl filter_contam_full_from_mothur_groups.plx file.groups contam.full

use warnings;
use strict;

#declare var
my $i=0;
my $contam;
my $line;
my $id;

#declare array
my @groups;
my @contam;
my @line;

#declare hash
my %contam;

open (IN, "<", $ARGV[0]) || die "Error cannot open groups file: $!\n";
@groups = <IN>;
close IN;

open (CONTAM, "<", $ARGV[1]) || die "Error cannot open contam.full file: $!\n";
@contam = <CONTAM>;
close CONTAM;

#add contam to hash for easy checking
while ($contam[$i]) {
	$contam = $contam[$i];
	chomp $contam;

	$contam{$contam} = 1;

	$i++;
	$contam=();
}
$i=0;

open ( OUT, ">>", "groups2.txt") || die "Error cannot open outfile: $!\n";

while ($groups[$i]) {
	$line = $groups[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$id = $line[0];

	if (exists $contam{$id}) {
		$i++;
		next;
	}
	else {
		print OUT "$line\n";
	}
	$i++;
}
$i=0;
close OUT;

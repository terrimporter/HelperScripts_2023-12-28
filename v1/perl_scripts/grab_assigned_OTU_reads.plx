#!/usr/bin/perl
#Dec.15, 2010 by Terri Porter
#Script to take a list of IDs from MEGAN (select leaf->export->reads) and the original clustering file from Uclust results.fasta and grab all component reads for each OTU ID
#usage perl grab_assigned_OTU_reads.plx MEGAN.fasta Uclust.results.fasta

use strict;
use warnings;

#var
my $line;
my $id;
my $i=0;
my $id_to_match;
my $j=0;
my $flag=0;

#array
my @uclust;
my @ids;

open (MEGAN, "<", $ARGV[0]) || die ("Error, cannot open megan file to read: $!\n");

open (UCLUST,"<", $ARGV[1]) || die ("Error, cannot open uclust file to read: $!\n");
@uclust = <UCLUST>;
close UCLUST;

while (<MEGAN>) {
	$line = $_;
	chomp $line;
	if (/>/) {
		$line =~ />(\w{14})/;
		$id = $1;
		push (@ids, $id);
	}
}
close MEGAN;

#print "@ids\n";#test

open (OUT, ">>", "assigned_OTU_reads.txt") || die ("Error, cannot open outfile to write: $!\n");

while ($ids[$i]) {
	$id_to_match = $ids[$i];

	while ($uclust[$j]) {
		$line = $uclust[$j];
		chomp $line;

		if ($flag==0) {
			if ($line =~ /\*/) {
				if ($line =~ /$id_to_match/) {
					print OUT $line."\n";
					$flag=1;
				}
			}
		}
		elsif ($flag==1) {
			if ($line =~ /\*/) {
				$flag=0;
				if ($line =~ /$id_to_match/) {
					print OUT $line."\n";
					$flag=1;
				}
			}
			else {
				print OUT $line."\n";
			}
		}
		$j++;
	}
	$j=0;
	$i++;

}

close OUT;


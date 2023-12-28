#!/usr/bin/perl
#Terri Porter, Aug.27, 2010
#Script to grab ID and cluster size from uc dereplication file
#usage $perl parse_uc_cluster_size.plx < results.uc > uc.clusters

use strict;
use warnings;

#declare variables
my $line;
my $id_line;
my $cluster_size;
my $id;

#declare array
my @line;

print "ID\tCluster_size\n";

while (<>) {
	$line = $_;
	chomp $line;
	if ($line =~ /^C/) {
		@line = split (/\t/, $line);
		$id_line = $line[8];
		$cluster_size = $line[2];
			if ($id_line =~ /\w{14,16}\s+/) {
				$id_line =~ /(\w{14,16})/;
				$id = $1;
				print "$id\t$cluster_size\n";
			}
	}
	else {
		next;
	}
}

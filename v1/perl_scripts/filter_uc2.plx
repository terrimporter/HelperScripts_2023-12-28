#!/usr/bin/perl
#Terri Porter, Aug. 26, 2010
#Script to filter sequences by cluster size
#usage $perl filter_uc.plx file.uc sorted_removenewline.fasta

use strict;
use warnings;

#declare variables
my $line;
my $cluster_size;
my $id_line;
my $id;
my $num;#test
#my $header;
my $i=0; #flag
my $j=0; #counter

#declare arrays
my @line;
my @ids;

filter();
#print "@ids\n"; #test
#print scalar(@ids); #test
print_filtered_fastas();

#subroutine to grab IDs of clusters with 3 or more reads

sub filter {
	open (IN,'<', $ARGV[0]) || die ("Error: $!\n");

	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^C/) {
			@line = split (/\t/, $line);
			$cluster_size = $line[2];
			$id_line = $line[8];
			if ($cluster_size >= 3) {  ###EDIT READ FREQUENCY FILTER HERE###
				$id_line =~ /\|(\w{14})\s+/;
				$id = $1;
				push (@ids, $id);
			}
		}
	}
close IN;
}

#subroutine to search fasta file according to a filtered list of ids from above

sub print_filtered_fastas {
	open (IN2, '<', $ARGV[1]) || die ("Error: $!\n");
	open (OUT, '>>', 'filtered_fasta.txt') || die ("Error: $!\n");
#print "@ids\t"; #test		
	while (<IN2>) {
		$line = $_;
		chomp $line;
		if ($ids[$j]) {
			if ($i == 0) {
				$id = $ids[$j];
				#print "$id\tTEST\t"; #test
				if ($line =~ /$id/) {
					print OUT $line."\n";
					$i = 1;
					next;
				}
				else {
					next;
					$i=0;
				}
			}
			elsif ($i == 1) {
				print OUT $line."\n";
				$i=0;
			}
			$j++;
		}
	}
close IN2;
close OUT;
}

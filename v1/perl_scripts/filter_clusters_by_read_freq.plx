#!/usr/bin/perl
#July 13, 2011 by Terri Porter
#Script to filter usearch .uc file by read frequency, then run blast
#do this after manually running FungalITSextractor, followed by parsing using sort_extracted_by_qual.plx, then manually running usearch --cluster file.fasta --uc results.uc --id 0.97 --nofastalign --rev --usersort --seedsout seed.fasta
#usage filter_clusters_by_read_freq.plx 97.uc file.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $cluster_size;
my $minClusterSize = 2; ####Edit minimum cluster size allowed here ####
my $id_line;
my $id;
my $j=0;
my $k=0;
my $l;
my $m=0;
my $header;
my $sequence;

#declare array
my @in;
my @line;
my @ids;
my @in2;
my @filtered;

open (IN,"<",$ARGV[0]) || die ("Error cannot read from results uc file: $!\n");
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
        chomp $line;
        
	if ($line =~ /^C/) {
        	@line = split (/\t/,$line);
                $cluster_size = $line[2];
                $id_line = $line[8];
                if ($cluster_size >= $minClusterSize) {
                	$id_line =~ /^(\w{14,16})/;
                        $id = $1;
                        push (@ids, $id);
                }
        }
	$i++;
}

open (IN2,"<",$ARGV[1]) || die ("Error cannot read from sorted fna file: $!\n");
@in2 = <IN2>;
close IN2;

while ($in2[$j]) {
	$line = $in2[$j];
        chomp $line;

	while ($ids[$k]) {
		$id = $ids[$k];

		if ($line =~ /^>/){
			$header = $line;
        		if ($header =~ /$id/) {
                		push(@filtered, $header);
				$l = $j+1;
				$sequence = $in2[$l];
				chomp $sequence;
				push(@filtered, $sequence);
			}
        	}
		$k++;
	}
	$k=0;
	$j++;
}

open (OUT, ">>", "sorted.fna.filtered") || die ("Error cannot write to sorted fna filtered file: $!\n");

while ($filtered[$m]) {
	$line = $filtered[$m];
	chomp $line;
	print OUT "$line\n";
	$m++;
}
close OUT;

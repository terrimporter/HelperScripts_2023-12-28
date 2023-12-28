#!/usr/bin/perl
#April 12, 2013 by Terri Porter
#Script to run Usearch6.0.307 OTU clustering for Nagissa's ITS seqs
#Both.fasta needs to be in current directory
#usearch6.0.307_i86linux32 needs to be in PATH
#usage perl run_usearch_OTU_clustering.plx

use strict;
use warnings;

#declare var
my $result;

#declare array
my @output;

$result = system('usearch6.0.307_i86linux32 -derep_fulllength Both.fasta -sizeout -output derep.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem derep.fa -id 0.99 -centroids denoised.fa -sizein -sizeout');

$result = system('usearch6.0.307_i86linux32 -sortbysize denoised.fa -output sorted.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem sorted.fa -id 0.97 -sizein -sizeout -centroids otus.fa');

$result = system('usearch6.0.307_i86linux32 -sortbysize otus.fa -minsize 2 -output otus_minsize2.fa');

#count number of seqs in each file and print to screen for easy recording

$result = qx(grep ">" Both.fasta | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" derep.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" denoised.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" sorted.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" otus.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" otus_minsize2.fa | wc -l);
chomp $result;
print $result."\t";

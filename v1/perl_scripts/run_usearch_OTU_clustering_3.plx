#!/usr/bin/perl
#May 14, 2013 edited to work with GC_40plus.fasta file instead
#May 8, 2013 edited to run clustering for Tara's baits
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

#$result = system('usearch6.0.307_i86linux32 -derep_fulllength Both.fasta -sizeout -output derep.fa');

#$result = system('usearch6.0.307_i86linux32 -cluster_smallmem derep.fa -id 0.99 -centroids denoised.fa -sizein -sizeout');

#$result = system('usearch6.0.307_i86linux32 -sortbysize denoised.fa -output sorted.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.99 -sizeout -centroids 99.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.98 -sizeout -centroids 98.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.97 -sizeout -centroids 97.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.96 -sizeout -centroids 96.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.95 -sizeout -centroids 95.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.94 -sizeout -centroids 94.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.93 -sizeout -centroids 93.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.92 -sizeout -centroids 92.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.91 -sizeout -centroids 91.fa');
$result = system('usearch6.0.307_i86linux32 -cluster_smallmem GC_40plus.fasta -id 0.90 -sizeout -centroids 90.fa');

#$result = system('usearch6.0.307_i86linux32 -sortbysize otus.fa -minsize 2 -output otus_minsize2.fa');

#count number of seqs in each file and print to screen for easy recording

$result = qx(grep ">" 99.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 98.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 97.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 96.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 95.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 94.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 93.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 92.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 91.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" 90.fa | wc -l);
chomp $result;
print $result."\n";

#$result = qx(grep ">" derep.fa | wc -l);
#chomp $result;
#print $result."\t";

#$result = qx(grep ">" denoised.fa | wc -l);
#chomp $result;
#print $result."\t";

#$result = qx(grep ">" sorted.fa | wc -l);
#chomp $result;
#print $result."\t";

#$result = qx(grep ">" otus.fa | wc -l);
#chomp $result;
#print $result."\t";

#$result = qx(grep ">" otus_minsize2.fa | wc -l);
#chomp $result;
#print $result."\t";

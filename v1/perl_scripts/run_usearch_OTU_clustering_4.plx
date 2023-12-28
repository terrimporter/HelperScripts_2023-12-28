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

#process ITS1 locus first
$result = system('cat ITS1_F.fasta ITS1_R.fasta > ITS1.fasta');

$result = system('usearch6.0.307_i86linux32 -derep_fulllength ITS1.fasta -sizeout -output ITS1_derep.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem ITS1_derep.fa -id 0.99 -centroids ITS1_denoised.fa -sizein -sizeout');

$result = system('usearch6.0.307_i86linux32 -sortbysize ITS1_denoised.fa -output ITS1_denoised.fas');

$result = system('usearch6.0.307_i86linux32 -uchime_denovo ITS1_denoised.fas -nonchimeras ITS1_nonch_denovo.fa -chimeras ITS1_ch_denovo.fa');

$result = system('usearch6.0.307_i86linux32 -uchime_ref ITS1_denoised.fas -db ITS1_nonch_denovo.fa -strand plus -nonchimeras ITS1_nonch_ref.fa -chimeras ITS1_ch_ref.fa');

$result = system('usearch6.0.307_i86linux32 -sortbysize ITS1_nonch_ref.fa -output ITS1_sorted.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem ITS1_sorted.fa -id 0.99 -sizein -sizeout -centroids ITS1_99otus.fa');

$result = system('usearch6.0.307_i86linux32 -sortbysize ITS1_99otus.fa -minsize 2 -output ITS1_99otus_minsize2.fa');

#now process ITS2 locus
$result = system('cat ITS2_F.fasta ITS2_R.fasta > ITS2.fasta');

$result = system('usearch6.0.307_i86linux32 -derep_fulllength ITS2.fasta -sizeout -output ITS2_derep.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem ITS2_derep.fa -id 0.99 -centroids ITS2_denoised.fa -sizein -sizeout');

$result = system('usearch6.0.307_i86linux32 -sortbysize ITS2_denoised.fa -output ITS2_denoised.fas');

$result = system('usearch6.0.307_i86linux32 -uchime_denovo ITS2_denoised.fas -nonchimeras ITS2_nonch_denovo.fa -chimeras ITS2_ch_denovo.fa');

$result = system('usearch6.0.307_i86linux32 -uchime_ref ITS2_denoised.fas -db ITS2_nonch_denovo.fa -strand plus -nonchimeras ITS2_nonch_ref.fa -chimeras ITS2_ch_ref.fa');

$result = system('usearch6.0.307_i86linux32 -sortbysize ITS2_nonch_ref.fa -output ITS2_sorted.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem ITS2_sorted.fa -id 0.99 -sizein -sizeout -centroids ITS2_99otus.fa');

$result = system('usearch6.0.307_i86linux32 -sortbysize ITS2_99otus.fa -minsize 2 -output ITS2_99otus_minsize2.fa');

#count number of seqs in each file and print to screen for easy recording
$result = qx(grep ">" ITS1.fasta | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_derep.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_denoised.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_nonch_denovo.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_ch_denovo.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_nonch_ref.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_ch_ref.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_99otus.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS1_99otus_minsize2.fa | wc -l);
chomp $result;
print $result."\n";

$result = qx(grep ">" ITS2.fasta | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_derep.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_denoised.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_nonch_denovo.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_ch_denovo.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_nonch_ref.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_ch_ref.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_99otus.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" ITS2_99otus_minsize2.fa | wc -l);
chomp $result;
print $result."\n";


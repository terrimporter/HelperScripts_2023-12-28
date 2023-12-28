#!/bin/bash
#Oct. 24, 2012 edited to run many small blast jobs on info nodes
#Oct. 18, 2012 by Terri Porter
#Script to automatically run a bunch of files with my_qsub

DIR="/home/terri/yersinia/R1_trimmed_MIN_MAX/all_blastn"
for i in $(ls $DIR | grep ".out");
do
	my_qsub -l "info101,info102,info103,info104,info105,info108,info109,info110,info111,info112" "blastn -query $DIR/$i -out $DIR/$i.blastn -outfmt 0 -task megablast -db refseq_genomic -word_size 28 -evalue '1e-10' -num_descriptions 100 -num_alignments 100"
done

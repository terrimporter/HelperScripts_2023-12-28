#!/bin/zsh
#Script to grab all paired end primer trimmed files, loop through, run seqprep multiple to pair the reads
#USAGE zsh runcutadapt.sh


for r1 in *R1.fastq.V3trimmed.fastq
do

echo $r1

base=$r1[1,4]
suffix=_R2.fastq.V4trimmed.fastq

r2=$base$suffix

echo $r2

#Pair 16S v3-v4 16S_V3_F and 16S_V4_R primer trimmd reads
seqprep -f $r1 -r $r2 -1 ${r1}.out -2 ${r2}.out -q 20 -s ${base}.16Spaired.fastq.gz -E ${base}.16Spaired.aln.gz -o 25

done
	


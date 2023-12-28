#!/bin/zsh
#Script to grab all paired end primer trimmed files, loop through, run seqprep multiple to pair the reads
#USAGE zsh runseqprep_BR5.sh


for r1 in *R1.fastq
do

echo $r1

base=$r1[1,4]
suffix=_R2.fastq

r2=$base$suffix

echo $r2

#Pair CO1 BR5 raw reads
seqprep -f $r1 -r $r2 -1 ${r1}.out -2 ${r2}.out -q 20 -s ${base}.paired.fastq.gz -E ${base}.paired.aln.gz -o 25

done
	


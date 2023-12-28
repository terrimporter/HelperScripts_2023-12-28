#!/bin/zsh
#Script to grab all paired end primer trimmed files, loop through, run seqprep multiple to pair the reads
#USAGE zsh runcutadapt.sh


for r1 in *R1.fastq.FFtrimmed.fastq
do

echo $r1

base=$r1[1,4]
suffix=_R2.fastq.FRtrimmed.fastq

r2=$base$suffix

echo $r2

#Pair Fish Mini_SH-F fish_miniF_F and fish_miniF_R primer trimmd reads
seqprep -f $r1 -r $r2 -1 ${r1}.out -2 ${r2}.out -q 20 -s ${base}.miniFpaired.fastq.gz -E ${base}.miniFpaired.aln.gz -o 25

done
	


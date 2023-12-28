#!/bin/zsh
#Script to grab all paired end primer trimmed files, loop through, run seqprep multiple to pair the reads
#USAGE zsh runcutadapt.sh


for r1 in *R1.fastq.EFtrimmed.fastq
do

echo $r1

base=$r1[1,4]
suffix=_R2.fastq.ERtrimmed.fastq

r2=$base$suffix

echo $r2

#Pair Fish Mini_SH-E fish_miniE_F and fish_miniE_R primer trimmd reads
seqprep -f $r1 -r $r2 -1 ${r1}.out -2 ${r2}.out -q 20 -s ${base}.miniEpaired.fastq.gz -E ${base}.miniEpaired.aln.gz -o 25

done
	


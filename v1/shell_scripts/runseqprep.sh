#!/bin/zsh
#Oct.3, 2016 by Terri Porter
#Script to grab all paired end seq files, loop through, run seqprep multiple to pair the reads
#USAGE zsh runseqprep.sh

NR_CPUS=10
count=0

for r1 in *_R1.fastq
do

echo $r1

base=${r1%%_R1.fastq}
suffix=_R2.fastq

r2=$base$suffix

echo $r2

#Pair raw reads
seqprep -f $r1 -r $r2 -1 ${r1}.out -2 ${r2}.out -q 20 -s ${base}.paired.fastq.gz -o 25 &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"

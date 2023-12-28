#!/bin/zsh
#Script to count number of lines in a bunch of files, copy to excel to count number of fastq seqs in a file
#USAGE zsh countlines.sh


for r1 in *.fastq
do

#echo $r1

base=$r1[1,3]
#suffix=_R2.fastq.Etrimmed

#r2=$base$suffix

#echo $r2

#wc -l trimmed files for CO1 BE amplicon, B, E
numlinesr1=$(wc -l $r1 | awk '{print $1}')
#numlinesr2=$(wc -l $r2 | awk '{print $1}')

echo $base'\t'$numlinesr1

done
	


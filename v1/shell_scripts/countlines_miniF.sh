#!/bin/zsh
#Script to count number of lines in a bunch of files, copy to excel to count number of fastq seqs in a file
#USAGE zsh countlines.sh


for r1 in *R1.fastq.FFtrimmed
do

#echo $r1

base=$r1[1,4]
suffix=_R2.fastq.FRtrimmed

r2=$base$suffix

#echo $r2

#wc -l trimmed files for Fish Mini_SH-F amplicon, Fish_miniF_F, Fish_miniF_R
numlinesr1=$(wc -l $r1 | awk '{print $1}')
numlinesr2=$(wc -l $r2 | awk '{print $1}')

echo $base'\t'$numlinesr1'\t'$numlinesr2

done
	


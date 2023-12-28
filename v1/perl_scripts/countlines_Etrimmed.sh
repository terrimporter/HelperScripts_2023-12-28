#!/bin/zsh
#Script to count number of lines in a bunch of files, copy to excel to count number of fastq seqs in a file
#USAGE zsh countlines.sh


for f in *.Etrimmed
do

base=$f[1,5]

#zcat | wc -l paired files
numlines=$(wc -l $f | awk '{print $1}')

echo $base'\t'$numlines

done
	


#!/bin/zsh
#Script to count number of fasta headers in a bunch of files, copy to excel
#USAGE zsh countfasta_uniques.sh


for f in *.uniques
do

base=$f[1,4]

#grep ">" filename | wc -l | awk '{print $1}'
numlines=$(grep ">" $f | wc -l | awk '{print $1}')

echo $base'\t'$numlines

done
	


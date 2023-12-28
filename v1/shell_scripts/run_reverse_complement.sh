#!/bin/zsh
#Script to run reverse_complement.plx on a directory of fasta files
#USAGE zsh run_reverse_complement.sh

for f in *.centroids
do

base=$f[1,6]

perl reverse_complement.plx $f

echo $base

done
	


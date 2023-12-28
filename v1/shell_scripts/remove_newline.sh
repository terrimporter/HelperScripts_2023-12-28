#!/bin/zsh
#Oct. 13, 2016 by Terri Porter
#Script to run remove_newline.plx on a directory of .centroids3 files
#USAGE zsh rename_fasta.sh

for f in *.centroids3
do

base=$f[1,5]

perl remove_newline.plx < $f > $f.fasta

echo $base

done
	


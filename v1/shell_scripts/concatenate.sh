#!/bin/zsh
#Oct. 13, 2016 by Terri Porter
#Script to concatenate a directory full of .blasn files so that they can be imported into MEGAN in one step
#Note that sample names need to be incorporated into fasta headers BEFORE blast to allow them to be disentangled later
#USAGE zsh concatenate.sh

for f in *.blastn
do

base=$f[1,5]

cat $f >> cat.blastn

echo $base

done
	


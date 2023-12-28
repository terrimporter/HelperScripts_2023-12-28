#!/bin/zsh
#Script to add sample name to fasta header
#USAGE zsh rename_fasta.sh

#echo sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

for f in *.centroids3
do

base=$f[1,5]

sed -i -- 's,>,>'"$base"'_,' "$f"

echo $base

done
	


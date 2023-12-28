#!/bin/zsh
#Script to add double quotes around each label in labels.txt for usearch
#USAGE zsh addquotes.sh

for f in *.txt
do

base=$f[1,6]

#perl fasta_stats_zsh.plx infile
sed 's/^/"/g; s/$/"/g' $f > $f.sed

echo $base

done
	


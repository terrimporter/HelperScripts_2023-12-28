#!/bin/zsh
#Script to run fasta_stats_zsh.plx on a directory of fasta files
#USAGE zsh run_fastastats.sh

echo sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

for f in *.centroids3
do

base=$f[1,4]

#perl fasta_stats_zsh.plx infile
stats=$( perl fasta_stats_zsh.plx $f )

echo $base'\t'$stats

done
	


#!/bin/zsh
#Script to run fasta_stats_zsh.plx on a directory of fasta files
#USAGE zsh run_fastastats.sh

echo primer'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

for f in *.fasta
do

#base=$f[1,4]

#echo $f

base=$(echo $f | cut -d"." -f1)

#perl fasta_stats_zsh.plx infile
stats=$( perl fasta_stats_zsh.plx $f )

echo $base'\t'$stats

done
	


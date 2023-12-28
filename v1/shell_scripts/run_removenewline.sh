#!/bin/zsh
#Script to run fasta_stats_zsh.plx on a directory of fasta files
#USAGE zsh run_fastastats.sh

#echo sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

for f in *.centroids3
do

base=$f[1,4]

###EDIT FILLE NAME AS NEEDED HERE###
fasta=$base.BE.fasta

#perl remove_newline_zsh.plx infile outfile
stats=$( perl remove_newline_zsh.plx $f $fasta)

#echo $base'\t'$stats

done
	


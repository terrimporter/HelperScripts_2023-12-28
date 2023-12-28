#!/bin/zsh
#Script to run check_seq_length.plx on a directory of fasta files
#USAGE zsh run_check_seq_length.sh

#echo primer'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

for f in *.fasta200
do

#base=$f[1,4]

#echo $f

base=$(echo $f | cut -d"." -f1)

#perl fasta_stats_zsh.plx infile
perl check_seq_length.plx $f

echo $base'\n'

done
	


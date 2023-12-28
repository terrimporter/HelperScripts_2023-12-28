#!/bin/zsh
#Jan. 30, 2017 edit to run jobs in parallel using zsh
#Script to run fasta_stats_zsh.plx on a directory of fasta files
#USAGE zsh run_fastastats.sh

echo sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

NR_CPUS=2
count=0

for f in *.uniques
do

stats $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"

#!/bin/zsh
#Jan. 30, 2017 edit to run jobs in parallel using zsh
#Script to run filter_by_orf.plx on a directory of fasta files
#USAGE zsh run_filter_by_orf.sh

NR_CPUS=10
count=0

for f in *.fasta
do

filterfasta_by_vertorf $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"

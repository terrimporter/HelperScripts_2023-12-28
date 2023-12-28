#!/bin/zsh
#Mar. 29/17 run ORFfinder on directory of files
#Jan. 30, 2017 edit to run jobs in parallel using zsh
#Script to run fasta_stats_zsh.plx on a directory of fasta files
#USAGE zsh run_fastastats.sh

NR_CPUS=10
count=0

for f in *.fasta
do

ORFfinder -in $f -g 2 -ml 309 > $f.vert &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"

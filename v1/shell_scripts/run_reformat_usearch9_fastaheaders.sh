#!/bin/zsh
#Feb. 6, 2017 edit to run reformat_usearch9_fastaheaders.sh on a directory of files
#Jan. 30, 2017 edit to run jobs in parallel using zsh
#USAGE zsh run_reformat_usearch9_fastaheaders.sh

NR_CPUS=4
count=0

for f in *.denoised
do

perl reformat_usearch9_fastaheaders.plx $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"

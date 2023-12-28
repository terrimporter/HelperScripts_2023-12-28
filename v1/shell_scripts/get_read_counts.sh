#!/bin/zsh
#Jan. 23, 2017 spread jobs across multiple cores
#Script to get read counts from fasta headers for all *.centroids3 in a directory
#USAGE sh get_read_counts.sh

#ngs-workflow only has 8 cores so don't run more than 7 at a once
NR_CPUS=4
count=0

echo -e 'sample\treadcount'

for f in *.sort2
do

base=${f%%.fastq.Ftrimmed.fasta.uniques.sort2*}
#echo $base

readcount=$(grep ">" $f | awk 'BEGIN {FS=";"} {print $2}' | sed 's/size=//g' | awk '{sum+=$1} END {print sum}')
#echo $readcount

echo -e $base'\t'$readcount

let count+=1 
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"


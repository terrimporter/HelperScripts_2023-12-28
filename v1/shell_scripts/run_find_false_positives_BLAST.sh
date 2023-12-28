#!/bin/bash
#Jun. 30/17 for each primer, check every merged blast out file with qcov against reference 
#USAGE sh run_find_false_positives_BLAST.sh

#ngs-workflow only has 8 cores so don't run more than 7 at a once
NR_CPUS=7
count=0

echo -e Primer'\t'TotalNumQueries'\t'truepositives'\t'falsenegatives'\t'truenegatives'\t'falsepositives

for f in *_merged_qcov.txt
do

base=${f%%_merged_qcov.txt}
#echo $base

find_false_positives $f mytrainseq.fasta &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"


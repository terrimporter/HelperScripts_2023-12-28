#!/bin/bash
#Jan. 23, 2017 by Terri Porter
#Script to add sample names and amplicon to headers then concatenate for easy import and sorting in excel
#USAGE sh rename_concatenate.sh

for f in *.out
do

base=${f%%.paired.fastq.Ftrimmed.fastq.fasta.uniques.sort2.centroids3.out*}

sed -i "s/^/$base;/g" "$f"

cat $f >> cat.out

echo $base

done

echo 'All jobs are done.'

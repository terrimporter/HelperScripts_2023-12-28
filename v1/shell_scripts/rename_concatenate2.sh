#!/bin/bash
#Jan. 25, 2017 edit to work with fasta files instead of BLAST or RDP output
#Jan. 23, 2017 by Terri Porter
#Script to add sample names and amplicon to headers then concatenate for easy import and sorting in excel
#USAGE sh rename_concatenate.sh

for f in *.centroids3
do

base=${f%%.paired.fastq.Ftrimmed.fastq.fasta.uniques.sort2.centroids3*}

sed -i "s/^>/>$base;/g" "$f"

cat $f >> cat.centroids3

echo $base

done

echo 'All jobs are done.'

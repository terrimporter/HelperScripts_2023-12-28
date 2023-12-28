#!/bin/bash
#Apr. 10/17, edit to work with denoised RDP outfiles
#Jan. 25, 2017 edit to work with fasta files instead of BLAST or RDP output
#Jan. 23, 2017 by Terri Porter
#Script to add sample names and amplicon to headers then concatenate for easy import and sorting in excel
#USAGE sh rename_concatenate.sh

for f in *.out
do

base=${f%%.updatedsizefasta.out*}

sed -i "s/^/$base\|/g" "$f"

cat $f >> cat.denoised.out

echo $base

done

echo 'All jobs are done.'

#!/bin/bash
#July 9, 2012 by Terri Porter
#Shell script to automtacally run cutadapt on a bunch of fastq.gz files

DIR="$HOME/cholera/R2_raw"

files="$(ls $DIR | grep .gz)"

for i in $files
do
	cutadapt -b AGATCGGAAGAGC -o $i.trimmed.gz $i -e 0 -m 20 -q 20 -z &
	wait
done

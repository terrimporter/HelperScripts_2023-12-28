#!/user/bin/bash

# Script to rename a bunch of symlinks

for f in *.fq 
do
	f2=$(echo $f | cut -d. -f5)
	f3="${f2}.fastq"
	mv $f $f3
	gzip $f3
done

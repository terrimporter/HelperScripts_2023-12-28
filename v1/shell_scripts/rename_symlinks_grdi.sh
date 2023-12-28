#!/user/bin/bash

# Script to rename a bunch of symlinks

for f in *.gz 
do
	f2=$(echo $f | cut -d- -f6)
	f3=$(echo $f2 | cut -d_ -f1,2,4)
	prefix="TR_"
	f4="${prefix}${f3}"
	f5=$(echo $f4 | cut -d. -f1)
	f6="${f5}.fastq.gz"
	mv "$f" "$f6"
done

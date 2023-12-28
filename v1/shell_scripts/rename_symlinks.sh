#!/user/bin/bash

# Script to rename a bunch of symlinks

for f in *.gz 
do
	f2=$(echo $f | cut -d. -f5-7)
	mv $f $f2
done

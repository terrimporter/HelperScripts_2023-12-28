#!/bin/bash
#September 30, 2009 by Terri Porter
#Shell script to automatically run muscle on a bunch of fasta files in a directory

DIR="/home/terri/from_estortho_082312"
for i in $(ls $DIR);
do
	OUT=$i.out
	muscle -in $DIR/$i -out $DIR/$OUT
done


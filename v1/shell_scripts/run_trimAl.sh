#!/bin/bash
#August 21, 2012 edited to run trimAl v1.4.rev9 instead
#September 30, 2009 by Terri Porter
#Script to automatically run Gblocks on a file containing muscle.out protein alignments

DIR="/home/terri/from_muscle_082312"
for i in $(ls $DIR);
do
	trimal -in $DIR/$i -out $DIR/$i.gappyout -gappyout 
done


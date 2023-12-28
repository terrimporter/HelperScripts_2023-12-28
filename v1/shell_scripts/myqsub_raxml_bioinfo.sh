#!/bin/bash
#Nov. 5, 2012 edited to add sleep delay to avoid missing submitted jobs
#Oct. 24, 2012 edited to run raxmlHPC on bioinfo nodes
#Oct. 18, 2012 by Terri Porter
#Script to automatically run a bunch of files with my_qsub

DIR="/home/terri/chytrid/terri_test17/raxml_1000/test_700s"
for i in $(ls $DIR | grep 'BS*');
do
	my_qsub -q "raxmlHPC -f d -s $DIR/$i -m PROTGAMMAWAGF -P $DIR/LGmodel -n $i.out"
	sleep 10
done

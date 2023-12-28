#!/bin/bash
#Oct. 18, 2012 by Terri Porter
#Script to automatically run a bunch of files with my_qsub

DIR="/home/terri/yersinia/R1_trimmed_MIN_MAX/taxonomy_crawl_new_blastx"
for i in $(ls $DIR | grep PART);
do
	my_qsub -l "info101,info102,info103,info104,info105,info108,info109,info110,info111,info112" "perl $DIR/taxonomy_crawl_map2.plx $DIR/$i $DIR/taxid.parsed.uniq"
done

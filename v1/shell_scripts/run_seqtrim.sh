#!/bin/bash
#Nov. 18, 2011 by Terri Porter
#script to run seqtrim.pl on a bunch of fasta and qual files organized in subdirectories

DIR="$HOME/Kapuskasing_raw/3_barcode_trimmed"
declare -i maxjobs=20
declare -i jobsrunning=0

run_seqtrim()
{
	fastafile="$(ls "$DIR/$dir" | grep fasta)"
	pathtofasta="$DIR/$dir/$fastafile"
	qualfile="$(ls "$DIR/$dir" | grep qual)"
	pathtoqual="$DIR/$dir/$qualfile"
	OUT="$DIR/$dir/seqtrim"
	seqtrim.pl -f $pathtofasta -q $pathtoqual -v --sm=50 --QV=20 --QW=10 --arrange vnq --o $OUT
}

run_jobs()
{
	if [ "$jobsrunning" -lt "$maxjobs" ]
	then
		echo "$dir running"
		run_seqtrim &
		let "jobsrunning += 1"
	else
		jobsrunning=0
		wait && run_jobs
	fi
}

for dir in `ls "$DIR/"`
do
	if [ -d "$DIR/$dir" ]
	then
		run_jobs
	fi
done

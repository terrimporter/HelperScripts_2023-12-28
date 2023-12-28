#!/bin/bash
#Oct.5, 2011 by Terri Porter
#Shell script to automatically run prottest on a bunch of .phy files

declare -i maxjobs=25
DIR="$HOME/prottest/ProtTest2.4"
declare -i jobsrunning=0
#echo $maxjobs
#echo $jobsrunning

prottest () 
{
	        OUT=$i.out
		runProtTest -i $i -o $OUT -+I+G F 
		}

testmax () {
	if [ "$jobsrunning" -lt "$maxjobs" ]
	then
		echo "$i"
		prottest &
		let "jobsrunning += 1"
	else 
		jobsrunning=0
		wait && testmax
	fi

}

files="$(ls $DIR | grep .phy)"

for i in $files
do
	testmax
done

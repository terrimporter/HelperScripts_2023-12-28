#!/usr/bin/bash
#Dec. 8, 2011 by Terri Porter
#Script to run phylip using abunch of raxml bootstrap files from run_raxml.sh and make_command_files.plx

DIR="$HOME/prottest/ProtTest2.4/prottest_infiles/relaxed/done/bootstrap"

files="$(ls $DIR | grep 'command$')"

for i in $files
do
	j=${i%%.*} #get id from file extension
	consense < $i &
	wait
	mv "outtree" "$j.outtree"
	mv "outfile" "$j.outfile"
done

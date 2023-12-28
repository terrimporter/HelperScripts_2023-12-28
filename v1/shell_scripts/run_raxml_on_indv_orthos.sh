#!/bin/bash
#Sept. 10, 2012 edited to run RtREV and Dayhoff models properly
#Oct.6, 2011 edited to run raxml jobs on a directory of infiles
#Oct.5, 2011 by Terri Porter
#Shell script to automatically run prottest on a bunch of .phy files

declare -i maxjobs=8
DIR="$HOME/prottest/ProtTest2.4/terri_test15/prottest_infiles/relaxed"
declare -i jobsrunning=0

set_up_command ()
{
	j=${i%%.*} #syntax ${variable%%pattern} trim longest match from end; get id from filename only
	prottest=$(grep $j id_model.map | awk '{print $2}') #get model+params
	model=${prottest%%+*} #get just model
	
	I=""
	F=""
	P=""

	if [ "$model" = "LG" ] #set up special case to use LGmodel
	then
		model="WAG"
		P="-P LGmodel"
	fi

	if [ "$model" = "RtREV" ] #change to CAPS for raxml
	then
		model="RTREV"
	fi

	if [ "$model" = "Dayhoff" ] #change to CAPS for raxml
	then 
		model="DAYHOFF"
	fi
	
	if [ "$model" = "DCMut" ] #change to CAPS for raxml
	then 
		model="DCMUT"
	fi

	if [ "$model" = "Blosum62" ] #change to CAPS for raxml
	then 
		model="BLOSUM62"
	fi


	if [[ "$prottest" =~ '+G+F' ]]
	then
		F="F"
	fi

	if [[ "$prottest" =~ '+I+F' ]]
	then 
		I="I"
		F="F"
	fi

	if [[ "$prottest" =~ '+F$' ]]
	then 
		F="F"
	fi
	
#	if [[ "$prottest" =~ 'Blosum62+I' ]]
#	then
#		I="I"
#	fi
	
	#raxml amino acid models all include GAMMA whether I want it or not!!!

	raxmlHPC -f a -x 23456 -N 100 -s $i $P -m PROTGAMMA$I$model$F -n $j
}

raxml () 
{
	        #OUT=$i.out
		#runProtTest -i $i -o $OUT -+I+G F 
		set_up_command
		}

testmax () {
	if [ "$jobsrunning" -lt "$maxjobs" ]
	then
		#echo "filename $i"
		raxml &
		let "jobsrunning += 1"
	else 
		jobsrunning=0
		wait && testmax
	fi

}

files="$(ls $DIR | grep .relaxed)"

for i in $files
do
	testmax
done

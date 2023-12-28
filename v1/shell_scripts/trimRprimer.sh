#!/bin/zsh
#Nov.22, 2016 by Terri Porter edit to use a primer file
#Script to grab all paired end files, loop through, run cutadapt multiple times to grab sequences with certain sets of trimmed primers
#USAGE zsh trimFprimer.sh Fprimers.txt
#USAGE zsh trimRprimer.sh Rprimers.txt

fileName="$1"

while read -r line
do

primerName=$(echo $line | cut -f1)
primerSeq=$(echo $line | cut -f2)

echo "Processing" $primerName

#trim 5' primer
#cutadapt -g $primerSeq --discard-untrimmed testNBC.fasta > $primerName.fasta

#trim 3' primer
cutadapt -a $primerSeq --discard-untrimmed testNBC.fasta > $primerName.fasta


done <"$fileName"
	


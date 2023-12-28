#!/bin/zsh
#Script to grab all fasta files, loop through, run rdp classifier on each one
#USAGE zsh run_rdp_CO1v4_NCBI.sh

#only have ~50G RAM so don't run more than 5-6 process at once!!!
NR_CPUS=4
count=0

for f in *.denoised
do

echo $f

java -Xmx8g -jar /home/terri/rdp_classifier_2.12/dist/classifier.jar classify -t /home/terri/CO1Classifier/v4/mydata/mydata_trained/rRNAClassifier.properties -o ${f}.out $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	


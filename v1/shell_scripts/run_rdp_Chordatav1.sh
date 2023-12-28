#!/bin/zsh
#Script to grab all fasta files, loop through, run rdp classifier on each one
#ensure symlink to fasta files is in rdp/
#USAGE zsh run_rdp_Arthropodav2.sh

#only have ~50G RAM so don't run more than 5-6 process at once!!!
NR_CPUS=5
count=0

for f in *.fasta
do

echo $f

#Run cutadapt RDP classifier using COI Arthropoda v2 training set
java -Xmx8g -jar /home/terri/rdp_classifier_2.12/dist/classifier.jar classify -t /home/terri/ChordataClassifier/v1/mydata/mydata_trained/rRNAClassifier.properties -o ${f}.out $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	


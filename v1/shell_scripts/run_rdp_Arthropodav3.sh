#!/bin/zsh
#Script to grab all fasta files, loop through, run rdp classifier on each one
#ensure symlink to fasta files is in rdp/
#USAGE zsh run_rdp_Arthropodav3.sh

#only have ~50G RAM so don't run more than 5-6 process at once!!!
NR_CPUS=3
count=0

for f in BlackbrookF230fastas/*.fasta
do

echo $f

#Run cutadapt RDP classifier using COI Arthropoda v3 training set
java -Xmx8g -jar dist/classifier.jar classify -t /home/terri/ArthropodaClassifier/v3/mydata/mydata_trained/rRNAClassifier.properties -o ${f}.out $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	


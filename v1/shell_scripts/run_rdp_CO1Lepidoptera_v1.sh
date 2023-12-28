#!/bin/zsh
#Script to grab all fasta files, loop through, run rdp classifier on each one
#ensure symlink to fasta files is in rdp/
#USAGE zsh run_rdp_CO1Lepidoptera_v1.sh

#only have ~50G RAM so don't run more than 5-6 process at once!!!
NR_CPUS=4
count=0

for f in *.denoised2
do

echo $f

#Run cutadapt RDP classifier using COI v2 training set that includes Arthropoda v3, Chordata v2, CO1 v1 outgroup taxa
java -Xmx8g -jar /home/terri/rdp_classifier_2.12/dist/classifier.jar classify -t /home/terri/CO1Lepidoptera/v1/mydata/mydata_trained/rRNAClassifier.properties -o ${f}.out $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	


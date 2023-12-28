#!/bin/zsh
#Script to grab all fasta files, loop through, run rdp classifier on each one
#ensure symlink to fasta files is in rdp/
#USAGE zsh run_rdp_18Sv3_2.sh

#only have ~50G RAM so don't run more than 5-6 process at once!!!
NR_CPUS=4
count=0

for f in *.denoised2
do

echo $f

java -Xmx8g -jar /home/terri/rdp_classifier_2.12/dist/classifier.jar classify -t /home/terri/18S_Eukaryota/v3.2/mydata/mydata_trained/rRNAClassifier.properties -o ${f}.out $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	


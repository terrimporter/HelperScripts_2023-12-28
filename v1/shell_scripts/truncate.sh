#!/bin/zsh
#Nov. 22, 2016 by Terri Porter
#Script to run truncate_seqs2.plx on a directory of fasta files
#USAGE zsh truncate_seqs2.plx [3/5] [200]

for file in *.fasta
do

echo $file

TrimFromEnd="$1"
FragSize="$2"

perl truncate_seqs2.plx $file $TrimFromEnd $FragSize

done

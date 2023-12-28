#!/usr/bin/zsh

# Teresita M. Porter, June 3, 2019
# Script to check two directories of GenBank retrieved records and to keep the biggest GenBank one in a new dir
# Had to do this because retrieving records at different times sometimes recovers more data?!?!?!
# Usage zsh keep_biggest_genbank.sh dir1 dir2 dir3

echo "PWD: $PWD"

DIR1="$1"
DIR2="$2"
DIR3="$3"

cd $DIR1

# grab files in dir1
for f in *.gb
do
#echo " f: $f"

FILE1=$PWD"/"$f
#echo "FILE1: $FILE1"
SIZE1=$(stat --format=%s $FILE1)
#echo "SIZE1: $SIZE1"

cd ".."

# check with same file in dir2
FILE2=$PWD"/"$DIR2"/"$f
#echo "FILE2: $FILE2"

SIZE2=$(stat --format=%s $FILE2)
#echo "SIZE2: $SIZE2"

FILE3=$PWD"/"$DIR3"/"$f
#echo "FILE3: $FILE3"

if [ $SIZE1 -gt $SIZE2 ]
then
	cp -f "$FILE1" "$FILE3"
#	echo "first bigger"

elif [ $SIZE2 -gt $SIZE1 ]
then
	cp -f "$FILE2" "$FILE3"
#	echo "second bigger"

else
	cp -f "$FILE1" "$FILE3"
#	echo "same size"
fi

cd $DIR1

done

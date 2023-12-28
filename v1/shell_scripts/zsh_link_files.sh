#!/bin/zsh
#May 30, 2017 by Terri Porter
#script to create a directory of files using symbolic links only

echo "Please enter path containing original files (OMIT final '/'):"
read path

echo "Enter file extention to search for (ex. '.fastq.gz'):"
read ext

cwd=$(pwd)
echo $cwd

for f in $path/*$ext
        do
                base=${f##$path/}
				echo $path
				echo $base
				/bin/ln \-s $f $base
        done


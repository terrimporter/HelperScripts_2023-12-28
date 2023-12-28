#!/bin/bash
#May 30, 2017 by Terri Porter
#script to create a directory of files using symbolic links only

echo "Please enter path containing original files (OMIT final '/'):"
read path

cwd=$(pwd)
echo $cwd

for f in $path/*.fasta
        do
                base=${f//$path\//}
                ln -s $path/$base $cwd/$base
        done


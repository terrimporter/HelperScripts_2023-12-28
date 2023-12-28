# finds file sizes in current directory
du -xk | sort -n | tail -25 > du.txt

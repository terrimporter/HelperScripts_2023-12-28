# checks each user to see amount of disk space they are using
find . -type f -printf "%u  %s\n" \
  | awk '{user[$1]+=$2}; END{for(i in user) print i,user[i]}' > find.txt

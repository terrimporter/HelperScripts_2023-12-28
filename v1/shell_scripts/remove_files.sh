# Teresita M. Porter, July 19, 2021
# Script to move unneeded files into their own directory for now
# https://unix.stackexchange.com/questions/333458/shell-script-for-moving-selected-files-from-one-directory-to-another

prev_dir=/home/terri/Atlantic_2021/MetaWorks1.8.0/data
new_dir=/home/terri/Atlantic_2021/MetaWorks1.8.0/data_removed
cd $prev_dir
for i in `cat /home/terri/Atlantic_2021/MetaWorks1.8.0/filenames_to_remove.csv`
do
#   sed -i 's/\r$//' $i 
   echo $i
   cd $prev_dir
   mv $i $new_dir
done

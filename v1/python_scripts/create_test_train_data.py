# Teresita M. Porter, Feb. 13/23
# aim for 5-fold cross validation sets
# split Metazoa data into 5 sets
# test seqs should only contain Metazoa
# ref seqs should include eukaryote & prokaryote outgroups
# USAGE python3 create_test_train_data.py testNBC_2.fasta

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import sys, re, random

# file names provided as command line arguments
fasta_file = sys.argv[1] # input fasta file as arg 1

# create empty dictionaries
metazoa_dict = {}
outgroup_dict = {}

# create test sets from metazoa only
fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file),'fasta'))

for key in fasta_dict:
	if re.search("Metazoa", fasta_dict[key].description):
		metazoa_dict[key] = str(fasta_dict[key].seq)
	else:
		outgroup_dict[key] = str(fasta_dict[key].seq)

# size of each set for 5-fold cross validation
n = len(metazoa_dict)//5
#print(n)

# copy the dic for work
temp_dict = metazoa_dict.copy()
# first set of keys
keys1 = random.sample(list(temp_dict), n)
for key in keys1:
	del temp_dict[key]
# second set of keys
keys2 = random.sample(list(temp_dict),n)
for key in keys2:
	del temp_dict[key]
# third set of keys
keys3 = random.sample(list(temp_dict),n)
for key in keys3:
	del temp_dict[key]
# fourth set of keys
keys4 = random.sample(list(temp_dict),n)
for key in keys4:
	del temp_dict[key]
# fifth set of keys
keys5 = random.sample(list(temp_dict),n)
for key in keys5:
	del temp_dict[key]
# get outgroup keys
keysout = list(outgroup_dict.keys())

# create 5 new directories with the unique metazoan keys
dict1 = {key:fasta_dict[key] for key in keys1}
dict2 = {key:fasta_dict[key] for key in keys2}
dict3 = {key:fasta_dict[key] for key in keys3}
dict4 = {key:fasta_dict[key] for key in keys4}
dict5 = {key:fasta_dict[key] for key in keys5}
dictout = {key:fasta_dict[key] for key in keysout}

# 5 test files
with open("test1.fasta", "w") as handle:
	SeqIO.write(dict1.values(), handle, "fasta-2line")

with open("test2.fasta", "w") as handle:
	SeqIO.write(dict2.values(), handle, "fasta-2line")

with open("test3.fasta", "w") as handle:
	SeqIO.write(dict3.values(), handle, "fasta-2line")

with open("test4.fasta", "w") as handle:
	SeqIO.write(dict4.values(), handle, "fasta-2line")

with open("test5.fasta", "w") as handle:
	SeqIO.write(dict5.values(), handle, "fasta-2line")

# 5 training files
# merge the other dict together to create train sets
train1 = dict2 | dict3 | dict4 | dict5 | dictout
with open("train1.fasta", "w") as handle:
	SeqIO.write(train1.values(), handle, "fasta-2line")

train2 = dict1 | dict3 | dict4 | dict5 | dictout
with open("train2.fasta" ,"w") as handle:
	SeqIO.write(train2.values(), handle, "fasta-2line")

train3 = dict1 | dict2 | dict4 | dict5 | dictout
with open("train3.fasta", "w") as handle:
	SeqIO.write(train3.values(), handle, "fasta-2line")

train4 = dict1 | dict2 | dict3 | dict5 | dictout
with open("train4.fasta", "w") as handle:
	SeqIO.write(train4.values(), handle, "fasta-2line")

train5 = dict1 | dict2 | dict3 | dict4 | dictout
with open("train5.fasta", "w") as handle:
	SeqIO.write(train5.values(), handle, "fasta-2line")

# get outgroup keys
#keysout = list(outgroup_dict.keys())

#with open("outgroup.fasta", "w") as handle:
#	SeqIO.write(dictout.values(), handle, "fasta-2line")



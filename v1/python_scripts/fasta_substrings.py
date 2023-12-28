# Teresita M. Porter, March 7, 2023
# Script to subsample differenet lengths from a FASTA sequence
# USAGE python3 fasta_substrings.py test1.fasta 200

import sys
from Bio import SeqIO
from random import randint

# read in FASTA file
infile = sys.argv[1]

# read in substring length
substringLength = sys.argv[2]

# create outfile name
outfile = "test5_"+substringLength+".fasta"
#print(outfile)

# list of new records
seq_frags = []

# parse through original FASTA file
fasta_seqs = SeqIO.parse(open(infile),'fasta')
with open(outfile, 'w+') as out_file:
	# parse through each record
	for fasta in fasta_seqs:

		# sample a random sequence fragment from each record
		limit = len(fasta)
		start = randint(0,limit-int(substringLength))
		end = start+int(substringLength)
		frag = fasta.seq[start:end]
		# assign the new frag to the record
		fasta.seq = frag
		# add record to a new list
		seq_frags.append(fasta)
		
# write everything out to a new FASTA file
	SeqIO.write(seq_frags, out_file, "fasta-2line")

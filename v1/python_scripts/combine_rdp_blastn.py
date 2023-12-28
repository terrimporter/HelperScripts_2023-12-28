# Teresita M. Porter, March 14, 2023
# Script to combine rdp classifier results with top BLAST hit results
# USAGE python3 combine_rdp_blastn.py rdp.out top.blastn

import sys
import pandas as pd

rdpfile = sys.argv[1]
blastfile = sys.argv[2]

rdp_df = pd.read_csv(rdpfile, sep='\t', header=None)
rdp_df.columns = ['QueryID','Strand','Root','rRank','rBP',
	'Superkingdom','skRank','skBP',
	'Kingdom','kRank','kBP',
	'Phylum','pRank','pBP',
	'Class','cRank','cBP',
	'Order','oRank','oBP',
	'Family','fRank','fBP',
	'Genus','gRank','gBP',
	'Species','sRank','sBP']
#print(rdp_df.head())

blast_df = pd.read_csv(blastfile, sep='\t', header=None)
blast_df.columns = ['QueryID','SubjectID','SubjectScientificName',
	'PercentIdentity','AlignmentLength','ExpectValue','BitScore','QueryCoveragePerSubject']
#print(blast_df.head())

combined_df = rdp_df.merge(blast_df, on='QueryID')
#print(combined_df.head())

combined_df.to_csv('test.out.blastn', index=False)

# Teresita M. Porter, March 13/23
# need to calculate total assigned in each bin (not just correct)
# print final df directly to CSV, print status messages to stdout without buffering
# Add some status messages
# Fix problem saving whole results to file, don't truncate df to fit to screen, ugh
# Remove optional column 'presence' to shorten run time
# Script to parse classifier output to create table that lists number of correct assignment for each bootstrap bin
# USAGE python3 calc_table.py test1.fasta train1.fasta test1.out

import sys
import os
from Bio import SeqIO
import pandas as pd
import numpy as np

print("Grabbing names of infiles...", flush=True)
testfile = sys.argv[1]
trainfile = sys.argv[2]
classified = sys.argv[3]
# automatically generate an outfile name
basename = os.path.splitext(testfile)[0]
outfile = basename+".table"

# read in testfile used as queries
# head -100 test1.fasta > test_data.fasta
# initialize mapping file
print("Parsing the test.fasta file...", flush=True)
test_df = pd.DataFrame()
with open(testfile, "r") as f:
	for record in SeqIO.parse(f,"fasta"):
		# add all headers to a list
		header = record.description
		df = pd.DataFrame([header], columns=['header'])
		df = df["header"].str.split(" ", expand = True)
		test_df = pd.concat([test_df,df])

# add column names
test_df.columns = ["Accession", "Lineage"]

# split lineage into fields
test_df2 = test_df.copy()
test_df2 = test_df2["Lineage"].str.split(";", expand = True)

# add column names
test_df2.columns = ["Root","Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"]

print("Merging dataframes...", flush=True)
# put it together
test_df3 = pd.concat([test_df, test_df2],axis = 1, join = 'outer', ignore_index=False, sort=False)

print("Creating species dictionary...", flush=True)
# Test species
test_species_list = set(test_df3['Species'])
# create species dict
test_species = test_df3[["Accession","Species"]]
test_species_dict = test_species.set_index("Accession")["Species"].to_dict()

print("Creating genus dictionary...", flush=True)
# Test genus
test_genus_list = set(test_df3['Genus'])
# create genus dict
test_genus = test_df3[["Accession","Genus"]]
test_genus_dict = test_genus.set_index("Accession")["Genus"].to_dict()

print("Creating family dictionary...", flush=True)
# Test family
test_family_list = set(test_df3['Family'])
# create family dict
test_family = test_df3[["Accession","Family"]]
test_family_dict = test_family.set_index("Accession")["Family"].to_dict()

print("Creating order dictionary...", flush=True)
# Test order
test_order_list = set(test_df3['Order'])
# create order dict
test_order = test_df3[["Accession","Order"]]
test_order_dict = test_order.set_index("Accession")["Order"].to_dict()

print("Creating class dictionary...", flush=True)
# Test class
test_class_list = set(test_df3['Class'])
# create class dict
test_class = test_df3[["Accession","Class"]]
test_class_dict = test_class.set_index("Accession")["Class"].to_dict()

print("Creating phylum dictionary...", flush=True)
# Test phylum
test_phylum_list = set(test_df3['Phylum'])
# create phylum dict
test_phylum = test_df3[["Accession","Phylum"]]
test_phylum_dict = test_phylum.set_index("Accession")["Phylum"].to_dict()

print("Creating kingdom dictionary...", flush=True)
# Test kingdom
test_kingdom_list = set(test_df3['Kingdom'])
# create kingdom dict
test_kingdom = test_df3[["Accession","Kingdom"]]
test_kingdom_dict = test_kingdom.set_index("Accession")["Kingdom"].to_dict()

print("Creating superkingdom dictionary...", flush=True)
# Test superkingdom
test_superkingdom_list = set(test_df3['Superkingdom'])
# create superkingdom dict
test_superkingdom = test_df3[["Accession","Superkingdom"]]
test_superkingdom_dict = test_superkingdom.set_index("Accession")["Superkingdom"].to_dict()







# read in trainfile used for training
# head -100 train1.fasta > train_data.fasta
# initialize mapping file
print("Parsing the train.fasta file...", flush=True)
train_df = pd.DataFrame()
with open(trainfile, "r") as f:
	for record in SeqIO.parse(f,"fasta"):
		# add all headers to a list
		header = record.description
		df = pd.DataFrame([header], columns=['header'])
		df = df["header"].str.split(" ", expand = True)
		train_df = pd.concat([train_df,df])

# add column names
train_df.columns = ["Accession", "Lineage"]

# split lineage into fields
train_df2 = train_df.copy()
train_df2 = train_df2["Lineage"].str.split(";", expand = True)

# add column names
train_df2.columns = ["Root","Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"]

print("Merging dataframes...", flush=True)
# put it together
train_df3 = pd.concat([train_df, train_df2],axis = 1, join = 'outer', ignore_index=False, sort=False)

# read in classifier results
# head -10 test1.out > test_data.out
print("Parsing the classifier results file...", flush=True)
df = pd.read_csv(classified, sep='\t', header=None)
#print(df.head())

# add column names
df.columns =['Accession', 'Strand', 'Root', 'rRank', 'RootBP',
'Superkingdom','skRank','skBP',
'Kingdom','kRank','kBP',
'Phylum','pRank','pBP',
'Class','cRank','cBP',
'Order','oRank','oBP',
'Family','fRank','fBP',
'Genus','gRank','gBP',
'Species','sRank','sBP']

# check Species assignments
print("Checking species assignments...", flush=True)
correct = []
total = []
support = []
for index,row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Species']
	acc = row['Accession']
	bp = row['sBP']
	support.append(bp)

	# get the original species from test_data.fasta
	target = test_species_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

species_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
species_df['rank'] = "Species"

# check Genus assignments
print("Checking genus assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Genus']
	acc = row['Accession']
	bp = row['gBP']
	support.append(bp)

	# get the original genus from test_data.fasta
	target = test_genus_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

genus_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct', 'total','support'])
genus_df['rank'] = "Genus"

# check Family assignments
print("Checking family assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Family']
	acc = row['Accession']
	bp = row['fBP']
	support.append(bp)

	# get the original family from test_data.fasta
	target = test_family_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

family_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
family_df['rank'] = "Family"

# check Order assignments
print("Checking order assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Order']
	acc = row['Accession']
	bp = row['oBP']
	support.append(bp)

	# get the original order from test_data.fasta
	target = test_order_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

order_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
order_df['rank'] = "Order"

# check Class assignments
print("Checking class assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Class']
	acc = row['Accession']
	bp = row['cBP']
	support.append(bp)

	# get the original class from test_data.fasta
	target = test_class_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

class_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
class_df['rank'] = "Class"

# check Phylum assignments
print("Checking phylum assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Phylum']
	acc = row['Accession']
	bp = row['pBP']
	support.append(bp)

	# get the original phylum from test_data.fasta
	target = test_phylum_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

phylum_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
phylum_df['rank'] = "Phylum"

# check Kingdom assignments
print("Checking kingdom assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Kingdom']
	acc = row['Accession']
	bp = row['kBP']
	support.append(bp)

	# get the original kingdom from test_data.fasta
	target = test_kingdom_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

kingdom_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
kingdom_df['rank'] = "Kingdom"

# check Superkingdom assignments
print("Checking superkingdom assignments...", flush=True)
correct = []
total = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Superkingdom']
	acc = row['Accession']
	bp = row['skBP']
	support.append(bp)

	# get the original superkingdom from test_data.fasta
	target = test_superkingdom_dict.get(acc)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
		total.append(1)
	else:
		correct.append(0)
		total.append(1)

superkingdom_df = pd.DataFrame(list(zip(correct, total, support)),
		columns=['correct','total','support'])
superkingdom_df['rank'] = "Superkingdom"

# put the df's together
print("Binning support values results...", flush=True)
cat_df = pd.concat([species_df, genus_df, family_df,
order_df, class_df, phylum_df,
kingdom_df, superkingdom_df])
# reindex
cat_df = cat_df.reset_index()

# copy to new column
cat_df['support_bins'] = cat_df['support']

# map support to bins
cat_df.loc[cat_df.support <= 0.09, 'support_bins'] = "0_9"
cat_df.loc[(cat_df.support >= 0.10) & (cat_df.support <= 0.19), 'support_bins'] = "10_19"
cat_df.loc[(cat_df.support >= 0.20) & (cat_df.support <=0.29), 'support_bins'] = "20_29"
cat_df.loc[(cat_df.support >= 0.30) & (cat_df.support <=0.39), 'support_bins'] = "30_39"
cat_df.loc[(cat_df.support >= 0.40) & (cat_df.support <=0.49), 'support_bins'] = "40_49"
cat_df.loc[(cat_df.support >= 0.50) & (cat_df.support <=0.59), 'support_bins'] = "50_59"
cat_df.loc[(cat_df.support >= 0.60) & (cat_df.support <=0.69), 'support_bins'] = "60_69"
cat_df.loc[(cat_df.support >= 0.70) & (cat_df.support <=0.79), 'support_bins'] = "70_79"
cat_df.loc[(cat_df.support >= 0.80) & (cat_df.support <=0.89), 'support_bins'] = "80_89"
cat_df.loc[(cat_df.support >= 0.90) & (cat_df.support <=0.94), 'support_bins'] = "90_94"
cat_df.loc[(cat_df.support >= 0.95) & (cat_df.support <=1.0), 'support_bins'] = "95_100"

# pivot table to summarize correct & total assignments at each rank for each support bin
print("Building final summary table...")

# create common index to accomodate for any missing support_bins
pt = pd.pivot_table(cat_df, values=['correct','total'], index=['rank'],
                    columns=['support_bins'], aggfunc=np.sum)

pt.columns = pt.columns.map('_'.join)

# get sums
pt['correct_Total'] = pt.filter(regex="correct_*").sum(axis=1)
pt['total_Total'] = pt.filter(regex="total_*").sum(axis=1)

# reorder rows
pt = pt.reindex(["Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])

# now fix up columns
# list of expected columns
colList = ['correct_0_9','total_0_9',
'correct_10_19','total_10_19',
'correct_20_29','total_20_29',
'correct_30_39','total_30_39',
'correct_40_49','total_40_49',
'correct_50_59','total_50_59',
'correct_60_69','total_60_69',
'correct_70_79','total_70_79',
'correct_80_89','total_80_89',
'correct_90_94','total_90_94',
'correct_95_100','total_95_100',
'correct_Total','total_Total']
# if each is present else create
ptFinal = pd.DataFrame()
for name in colList:
	if(name in pt.columns):
		ptFinal[name] = pt[name]
	else:
		ptFinal[name] = 0

# save as CSV
ptFinal.to_csv(outfile)

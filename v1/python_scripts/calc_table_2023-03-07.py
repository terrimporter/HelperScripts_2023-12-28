# Teresita M. Porter, March 7/23
# Remove optional column 'presence' to shorten run time
# Script to parse classifier output to create table that lists number of correct assignment for each bootstrap bin
# USAGE python3 calc_table.py test1.fasta train1.fasta test1.out > test1.table

import sys
from Bio import SeqIO
import pandas as pd
import numpy as np

testfile = sys.argv[1]
trainfile = sys.argv[2]
outfile = sys.argv[3]

# read in testfile used as queries
# head -100 test1.fasta > test_data.fasta
# initialize mapping file
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

# put it together
test_df3 = pd.concat([test_df, test_df2],axis = 1, join = 'outer', ignore_index=False, sort=False)

# Test species
test_species_list = set(test_df3['Species'])
# create species dict
test_species = test_df3[["Accession","Species"]]
test_species_dict = test_species.set_index("Accession")["Species"].to_dict()

# Test genus
test_genus_list = set(test_df3['Genus'])
# create genus dict
test_genus = test_df3[["Accession","Genus"]]
test_genus_dict = test_genus.set_index("Accession")["Genus"].to_dict()

# Test family
test_family_list = set(test_df3['Family'])
# create family dict
test_family = test_df3[["Accession","Family"]]
test_family_dict = test_family.set_index("Accession")["Family"].to_dict()

# Test order
test_order_list = set(test_df3['Order'])
# create order dict
test_order = test_df3[["Accession","Order"]]
test_order_dict = test_order.set_index("Accession")["Order"].to_dict()

# Test class
test_class_list = set(test_df3['Class'])
# create class dict
test_class = test_df3[["Accession","Class"]]
test_class_dict = test_class.set_index("Accession")["Class"].to_dict()

# Test phylum
test_phylum_list = set(test_df3['Phylum'])
# create phylum dict
test_phylum = test_df3[["Accession","Phylum"]]
test_phylum_dict = test_phylum.set_index("Accession")["Phylum"].to_dict()

# Test kingdom
test_kingdom_list = set(test_df3['Kingdom'])
# create kingdom dict
test_kingdom = test_df3[["Accession","Kingdom"]]
test_kingdom_dict = test_kingdom.set_index("Accession")["Kingdom"].to_dict()

# Test superkingdom
test_superkingdom_list = set(test_df3['Superkingdom'])
# create superkingdom dict
test_superkingdom = test_df3[["Accession","Superkingdom"]]
test_superkingdom_dict = test_superkingdom.set_index("Accession")["Superkingdom"].to_dict()








# read in trainfile used for training
# head -100 train1.fasta > train_data.fasta
# initialize mapping file
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

# put it together
train_df3 = pd.concat([train_df, train_df2],axis = 1, join = 'outer', ignore_index=False, sort=False)

# Trained species
#train_species_list = set(train_df3['Species'])
# Trained genus
#train_genus_list = set(train_df3['Genus'])
# Trained family
#train_family_list = set(train_df3['Family'])
# Trained order
#train_order_list = set(train_df3['Order'])
# Trained class
#train_class_list = set(train_df3['Class'])
# Trained phylum
#train_phylum_list = set(train_df3['Phylum'])
# Trained kingdom
#train_kingdom_list = set(train_df3['Kingdom'])
# Trained superkingdom
#train_superkingdom_list = set(train_df3['Superkingdom'])













# read in outfile, classifier results
# head -10 test1.out > test_data.out
df = pd.read_csv(outfile, sep='\t', header=None)
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
#presence = []
correct = []
support = []
for index,row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Species']
	acc = row['Accession']
	bp = row['sBP']
	support.append(bp)

	# get the original species from test_data.fasta
	target = test_species_dict.get(acc)
#	print("target: ",target,"asst: ",asst)

	# check if taxon present in reference set
#	if asst in train_species_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

species_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
species_df['rank'] = "Species"

#print(species_df)


# check Genus assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Genus']
	acc = row['Accession']
	bp = row['gBP']
	support.append(bp)

	# get the original genus from test_data.fasta
	target = test_genus_dict.get(acc)

	# check if taxon present in reference set
#	if asst in train_genus_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

genus_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
genus_df['rank'] = "Genus"

#print(genus_df)








# check Family assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Family']
	acc = row['Accession']
	bp = row['fBP']
	support.append(bp)

	# get the original family from test_data.fasta
	target = test_family_dict.get(acc)

	# check if taxon present in reference set
#	if asst in train_family_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

family_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
family_df['rank'] = "Family"

#print(family_df)







# check Order assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Order']
	acc = row['Accession']
	bp = row['oBP']
	support.append(bp)

	# get the original order from test_data.fasta
	target = test_order_dict.get(acc)

	# check if taxon present in reference set
#	if asst in train_order_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

order_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
order_df['rank'] = "Order"

#print(order_df)







# check Class assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Class']
	acc = row['Accession']
	bp = row['cBP']
	support.append(bp)

	# get the original class from test_data.fasta
	target = test_class_dict.get(acc)

	# check if present in reference set
#	if asst in train_class_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

class_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
class_df['rank'] = "Class"

#print(class_df)







# check Phylum assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Phylum']
	acc = row['Accession']
	bp = row['pBP']
	support.append(bp)

	# get the original phylum from test_data.fasta
	target = test_phylum_dict.get(acc)

	# check if taxon present in reference set
#	if asst in train_phylum_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

phylum_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
phylum_df['rank'] = "Phylum"

#print(phylum_df)









# check Kingdom assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Kingdom']
	acc = row['Accession']
	bp = row['kBP']
	support.append(bp)

	# get the original kingdom from test_data.fasta
	target = test_kingdom_dict.get(acc)

	# check if taxon present in reference set
#	if asst in train_kingdom_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

kingdom_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
kingdom_df['rank'] = "Kingdom"

#print(kingdom_df)









# check Superkingdom assignments
#presence = []
correct = []
support = []
for index, row in df.iterrows():
	# taxonomic assignment from test_data.out
	asst = row['Superkingdom']
	acc = row['Accession']
	bp = row['skBP']
	support.append(bp)

	# get the original superkingdom from test_data.fasta
	target = test_superkingdom_dict.get(acc)

	# check if taxon present in reference set
#	if asst in train_superkingdom_list:
#		presence.append(1)
#	else:
#		presence.append(0)

	# check if assignment is correct
	if asst == target:
		correct.append(1)
	else:
		correct.append(0)

superkingdom_df = pd.DataFrame(list(zip(correct, support)),
		columns=['correct','support'])
superkingdom_df['rank'] = "Superkingdom"

#print(superkingdom_df)








# put the df's together
cat_df = pd.concat([species_df, genus_df, family_df,
order_df, class_df, phylum_df,
kingdom_df, superkingdom_df])
# reindex
cat_df = cat_df.reset_index()
#print(cat_df.head(20))

# copy to new column
cat_df['support_bins'] = cat_df['support']
#with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#    print(cat_df)

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
#print(cat_df.head(20))

# pivot table to summarize correct assignments at each rank for each support bin
pt = pd.pivot_table(cat_df, values='correct', index=['rank'],
                    columns=['support_bins'], aggfunc=np.sum)
#print(pt)
#print(type(pt))




# count total number of correct assignments at each rank
pt.loc[:,'Correct'] = pt.sum(axis=1)
#print(pt)





# count total number of assignments at each rank
occur = cat_df.groupby(['rank']).size().to_dict()
occur_list = [occur['Species'], occur['Genus'], occur['Family'],
		occur['Order'], occur['Class'], occur['Phylum'],
		occur['Kingdom'], occur['Superkingdom']]
#print(occur_list)
occur_df = pd.DataFrame(occur_list, 
		index=['Species','Genus','Family','Order',
		'Class','Phylum','Kingdom','Superkingdom'])
occur_df.columns = ['Total']
occur_df.index.name = 'rank'

final_pt = pd.merge(pt, occur_df, on=['rank'])
final_pt = final_pt.loc[["Superkingdom", "Kingdom", "Phylum",
"Class","Order","Family","Genus","Species"],:]

# need to reset option_context so the output isn't truncated, ugh
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
	print(final_pt)


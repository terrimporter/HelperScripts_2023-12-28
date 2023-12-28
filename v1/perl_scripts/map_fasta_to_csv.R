# Teresita M. Porter, Sept. 12, 2019
# Map denoised ESV sequences to taxonomy table
# Usage Rscript map_fasta_to_csv.R cat.denoised taxonomy.csv

library(Biostrings)

# get infiles from command line
#args=(commandArgs(TRUE))

# read FASTA file
fasta <- readDNAStringSet("F230R.denoised")
seq <- as.data.frame(readDNAStringSet("F230R.denoised"))$x
fasta.df<-data.frame(header=names(fasta), seq=seq)

# read taxonomy matrix
tax<-read.table(file="GJ_COI_taxonomic_assignments_raw.csv", head=TRUE, sep=",")

# merge
merged<-merge(fasta.df,tax,by.x="header",by.y="Amplicon_GlobalESV")

# create lineage, finest level of well suppored resolution (99% correct at all ranks, except 95% at species rank)
# see https://github.com/terrimporter/CO1Classifier
merged$lineage <- paste(merged$Phylum, ";", merged$Class, ";", merged$Order, sep="")
merged$lineage <- ifelse(merged$fBP>=0.30, paste(merged$lineage, ";", merged$Family, sep=""), paste(merged$lineage, ";", sep=""))
merged$lineage <- ifelse(merged$gBP>=0.20, paste(merged$lineage, ";", merged$Genus, sep=""), paste(merged$lineage,";", sep=""))
merged$lineage <- ifelse(merged$sBP>=0.70, paste(merged$lineage, ";", merged$Species, sep=""), paste(merged$lineage, ";", sep=""))

# save
write.csv(merged, file = "F230R_matrix.csv", row.names=FALSE, quote=FALSE)


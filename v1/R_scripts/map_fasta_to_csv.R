# Teresita M. Porter, Sept. 12, 2019
# Map denoised ESV sequences to taxonomy table
# Usage Rscript map_fasta_to_csv.R cat.denoised taxonomy.csv

library(Biostrings)

# get infiles from command line
#args=(commandArgs(TRUE))

# read FASTA file
fasta <- readDNAStringSet("F230R.denoised")
seq <- as.data.frame(readDNAStringSet("F230R.denoised"))$x
fasta.merged<-data.frame(header=names(fasta), seq=seq)

# read taxonomy matrix
tax<-read.table(file="GJ_COI_taxonomic_assignments_raw.csv", head=TRUE, sep=",")

# merge
merged<-merge(fasta.merged,tax,by.x="header",by.y="Amplicon_GlobalESV")

# create lineage, finest level of well suppored resolution (99% correct at all ranks, except 95% at species rank)
# see https://github.com/terrimporter/CO1Classifier
merged$lineage <- paste(merged$Phylum, ";", merged$Class, ";", merged$Order, sep="")
merged$lineage <- ifelse(merged$fBP>=0.30, paste(merged$lineage, ";", merged$Family, sep=""), paste(merged$lineage, ";", sep=""))
merged$lineage <- ifelse(merged$gBP>=0.20, paste(merged$lineage, ";", merged$Genus, sep=""), paste(merged$lineage,";", sep=""))
merged$lineage <- ifelse(merged$sBP>=0.70, paste(merged$lineage, ";", merged$Species, sep=""), paste(merged$lineage, ";", sep=""))

merged[36:38] <- data.frame(do.call('rbind', strsplit(as.character(merged$SampleName),'_',fixed=TRUE)), stringsAsFactors = FALSE)
names(merged)[36:38] <- c("molecule","filter","PCRStep")

# Add bottle column based on filter (only interested in bottles 1-6)
merged$bottle <- ifelse(merged$filter %in%  c("1","9"), "1", NA)
merged$bottle <- ifelse(merged$filter %in%  c("2","10"), "2", merged$bottle)
merged$bottle <- ifelse(merged$filter %in%  c("3","11"), "3", merged$bottle)
merged$bottle <- ifelse(merged$filter %in%  c("4","12"), "4", merged$bottle)
merged$bottle <- ifelse(merged$filter %in%  c("5","13"), "5", merged$bottle)
merged$bottle <- ifelse(merged$filter %in%  c("6","14"), "6", merged$bottle)
merged$bottle <- ifelse(merged$filter %in%  c("7"), NA, merged$bottle)

# Add site column based on filter
merged$site <- ifelse(merged$filter %in%  c("1","9","2","10","3","11"), "A", "")
merged$site <- ifelse(merged$filter %in%  c("4","12","5","13","6","14"), "B", merged$site)

# save
write.csv(merged, file = "F230R_matrix.csv", row.names=FALSE, quote=FALSE)


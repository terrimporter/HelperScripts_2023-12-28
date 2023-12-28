library(vegan)
library(ggplot2)
library(ggpubr)
library(reshape)
library(splitstackshape)
library(scales)

# Create rarefaction curves for all phyla

## Read infile prepared by python script
A<-read.table(file="AllPhyla_LabCode_SiteCode_matrix.csv", head=TRUE, row.names=1, sep=",")
head(A)

## Only process labs sequenced at Guelph (drop Curtin), drop results for lab E
A2<-A[- grep("E_", rownames(A)),]
A3<-A2[- grep("G_", rownames(A2)),]
A4<-A3[- grep("H_", rownames(A3)),]
A5<-A4[- grep("I_", rownames(A4)),]
A6<-A5[- grep("J_", rownames(A5)),]
A7<-A6[- grep("K_", rownames(A6)),]
A8<-A7[- grep("L_", rownames(A7)),]

## Only process 4 sites (skip blanks and notemplate positives)
A9<-A8[- grep("_17", rownames(A8)),]
A10<-A9[- grep("_18", rownames(A9)),]
A11<-A10[- grep("_19", rownames(A10)),]
A12<-A11[- grep("_20", rownames(A11)),]

## Remove columns with only zeros
A_notnull<-A12[,colSums(A12) !=0]

## Remove rows with only zeros
A_notnull2<-A_notnull[rowSums(A_notnull) !=0,]

## Calculate 15th percentile for rrarefy function 
A_15percentile<-quantile(rowSums(A_notnull2), prob=0.15)

###################################################################
##### Plot rarefaction curves
###################################################################

# Print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile

pdf("RareCurves_AllPhyla.pdf")

A_rarecurveout<-rarecurve(A_notnull2, sample=A_15percentile, step=100, main="18S - All Eukaryota Phyla", xlab="Reads", ylab="Phyla", col="blue4", cex=0.6)

dev.off()

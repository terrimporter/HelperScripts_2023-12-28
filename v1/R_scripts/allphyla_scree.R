library(vegan)
library(ggplot2)
library(ggpubr)
library(reshape)
library(splitstackshape)
library(scales)

###################################################################
##### NEW FUNCTION FOR SCREE PLOT #################################
###################################################################

#new function edited from http://www.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm
nmds.scree<-function(x) {
	plot(rep(1,10),replicate(10,metaMDS(x, autotransform=F,k=1)$stress/100),xlim=c(1,nrow(x)),ylim=c(0,0.5),xlab="Number dimensions",ylab="Stress",main="NMDS stress plot")
	
	for(i in 1:(nrow(x)-2)) {
		points(rep(i+1,10),replicate(10,metaMDS(x, autotransform=F,k=i+1)$stress/100))
	}
}

#################

# Create stacked bar plots for target taxa, reads vs lab group, sites pooled

## Read infile prepared by python script
A<-read.table(file="AllPhyla_LabCode_SiteCode_matrix.csv", head=TRUE, row.names=1, sep=",")
head(A)
## Only process labs sequenced at Guelph (skip Curtin), throw out results for lab E
grepA<-A[grep("A_",rownames(A)),]
grepB<-A[grep("B_",rownames(A)),]
grepC<-A[grep("C_",rownames(A)),]
grepD<-A[grep("D_",rownames(A)),]
grepF<-A[grep("F_",rownames(A)),]
A2 <- rbind(grepA, grepB, grepC, grepD, grepF)

# Remove columns with only zeros
A_notnull<-A2[,colSums(A2) !=0]

# Remove rows with only zeros
A_notnull2<-A_notnull[rowSums(A_notnull) !=0,]

# Calculate 15th percentile for rrarefy function 
### Rarefy to 30th percentile here because too much data loss
A_15percentile<-quantile(rowSums(A_notnull2), prob=0.30)

###################################################################
##### Rarefy the dataset down to the 30th percentile
###################################################################

#Rarefy original matrix down to 30th percentile library size to normalize read depth across samples
A_df<-rrarefy(A_notnull2,sample=A_15percentile)

#Transform rarified matrix into a presence-absence matrix for Beta div calcns
A_df[A_df>0] <-1

###################################################################
##### Do scree plots
###################################################################

#Scree plots to determine number of dimensions to use for NMDS

pdf("allphyla_scree.pdf")
nmds.scree(A_df)
dev.off()

###################################################################
##### Create distance matrix for beta diversity analyses Phyla only
###################################################################

#sorensen dissimilarity based on presence/absence
sor<-vegdist(A_df, "bray", binary=TRUE)

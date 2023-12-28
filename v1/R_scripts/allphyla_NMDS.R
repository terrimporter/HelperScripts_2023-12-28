library(vegan)
require(graphics)
require(dendextend)
library(ggplot2)
library(ggpubr)
library(reshape)
library(splitstackshape)
library(scales)


set.seed(1234)

## Read infile prepared by python script
A<-read.table(file="AllPhyla_LabCode_SiteCode_matrix.csv", head=TRUE, row.names=1, sep=",")
#head(A)

## Only process labs sequenced at Guelph (skip Curtin), throw out results for lab E
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

# Remove columns with only zeros
A_notnull<-A12[,colSums(A12) !=0]

# Remove rows with only zeros
A_notnull2<-A_notnull[rowSums(A_notnull) !=0,]

# Calculate 15th percentile for rrarefy function 
### Rarefy to 30th percentile here because too much data loss
A_15percentile<-quantile(rowSums(A_notnull2), prob=0.30)

###################################################################
##### Rarefy the dataset down to the 15th percentile
###################################################################

#Rarefy original matrix down to 30th percentile library size to normalize read depth across samples
A_df<-rrarefy(A_notnull2,sample=A_15percentile)

#Transform rarified matrix into a presence-absence matrix for beta div calcns
A_df[A_df>0] <-1

###################################################################
##### Create distance matrix for beta diversity analyses Phyla only
###################################################################

#Do 2 dimensional NMDS, no environment file for now 
pdf("allphyla_NMDS_k2.pdf")
A_nmds2<-metaMDS(A_df, k=2,trymax=100)
print(row.names(A_df))

#color by sites (blue=Monterey, green=Victoria, orange=Hillary, purple=Waitemata)
cols=c(rep("blue4",4),rep("green4",4),rep("orange4",4),rep("purple4",4),
       rep("blue4",4),rep("green4",4),rep("orange4",4),rep("purple4",4),
       rep("blue4",4),rep("green4",3),rep("orange4",3),rep("purple4",4),
       rep("blue4",4),rep("green4",4),rep("orange4",4),rep("purple4",3),
       rep("blue4",3),rep("green4",4),rep("orange4",4),rep("purple4",4))

#shapes by labs (pch 15=A, 16=B, 17=C, 18=D, 19=F)
pch=c(rep(0,16),
      rep(1,16),
      rep(2,14),
      rep(3,15),
      rep(4,15))

#create dataframe to indicate groups
groups_df<-data.frame("names"=rownames(A_notnull2))
groups_df$lab_site <- groups_df$names
groups_df<-cSplit(groups_df,"lab_site",sep="_", direction="wide")
names(groups_df)[2]<-"lab"
names(groups_df)[3]<-"sitecode"

## Copy column
groups_df$Site<-groups_df$sitecode

## Rename numbered sites into their site names
groups_df$Site <- gsub('10', 'Hillary', groups_df$Site)
groups_df$Site <- gsub('11', 'Hillary', groups_df$Site)
groups_df$Site <- gsub('12', 'Hillary', groups_df$Site)
groups_df$Site <- gsub('13', 'Waitemata', groups_df$Site)
groups_df$Site <- gsub('14', 'Waitemata', groups_df$Site)
groups_df$Site <- gsub('15', 'Waitemata', groups_df$Site)
groups_df$Site <- gsub('16', 'Waitemata', groups_df$Site)
groups_df$Site <- gsub('1', 'Monterey', groups_df$Site)
groups_df$Site <- gsub('2', 'Monterey', groups_df$Site)
groups_df$Site <- gsub('3', 'Monterey', groups_df$Site)
groups_df$Site <- gsub('4', 'Monterey', groups_df$Site)
groups_df$Site <- gsub('5', 'Victoria', groups_df$Site)
groups_df$Site <- gsub('6', 'Victoria', groups_df$Site)
groups_df$Site <- gsub('7', 'Victoria', groups_df$Site)
groups_df$Site <- gsub('8', 'Victoria', groups_df$Site)
groups_df$Site <- gsub('9', 'Hillary', groups_df$Site)

#grab sites & scores
A_spp.sc <- scores(A_nmds2, display = "species")
A_site.sc <- scores(A_nmds2, display = "sites")

plot(A_nmds2, type="n", xlab="NMDS1", ylab="NMDS2", main="All Eukaryote Phyla (Guelph)", cex.axis=1.5, cex.lab=1.5, cex.main=2)
ordihull(A_nmds2,groups_df$Site, draw="polygon", border=NA, col=c("orange4","blue4","green4","purple4"))
points(A_site.sc[,1], A_site.sc[,2], col=cols, pch=pch)

legend("bottomleft",
        c("Monterey", "Victoria", "Hillary", "Waitemata", "A", "B", "C","D","F"),
        col=c("blue4","green4","orange4","purple4",rep("black",5)),
        pch=c(15,15,15,15,0,1,2,3,4,5),
        cex=1, pt.cex=1,
        bty="n")

dev.off()

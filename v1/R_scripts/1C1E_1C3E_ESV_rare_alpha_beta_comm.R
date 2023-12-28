library(vegan)
require(graphics)
require(dendextend)

###################################################################
##### EDIT VEGAN RARECURVE FUNCTION TO ADJUST XLAB AND YLAB #######
###################################################################

rarecurve_ESV<-function (x, step, sample,main="",xlab="Reads", ylab="ESVs",label=TRUE, col, lty, ...){
	x <- as.matrix(x)
	if (!identical(all.equal(x, round(x)), TRUE))
	stop("function accepts only integers (counts)")
	if (missing(col))
	col <- par("col")
	if (missing(lty))
	lty <- par("lty")
	tot <- rowSums(x)
	S <- specnumber(x)
	if (any(S <= 0)) {
		message("empty rows removed")
		x <- x[S > 0, , drop = FALSE]
		tot <- tot[S > 0]
		S <- S[S > 0]
	}
	nr <- nrow(x)
	col <- rep(col, length.out = nr)
	lty <- rep(lty, length.out = nr)
	out <- lapply(seq_len(nr), function(i) {
		n <- seq(1, tot[i], by = step)
		if (n[length(n)] != tot[i])
		n <- c(n, tot[i])
		drop(rarefy(x[i, ], n))
	})
	Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
	Smax <- sapply(out, max)
	plot(c(1, max(Nmax)), c(1, max(Smax)),main=main, xlab = xlab, ylab = ylab,type = "n", ...)
	if (!missing(sample)) {
		abline(v = sample)
		rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),y = z, xout = sample, rule = 1)$y)
		abline(h = rare, lwd = 0.5)
	}
	for (ln in seq_along(out)) {
		N <- attr(out[[ln]], "Subsample")
		lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)						    			   }
	if (label) {
		ordilabel(cbind(tot, S), labels = rownames(x), ...)
	}
	invisible(out)
}

###################################################################

#read infile prepared by python script
A<-read.csv(file="1C1E_1C3E_ESV_matrix.csv", head=TRUE, row.names=1)

#split out ILC and NZC datasets
ILC<-A[1:47,]
NZC<-A[48:94,]

#read infile metadata
ENV<-read.table(file="file.env", sep=",", head=TRUE, row.names=1)

#create metadata file for ILC and NZC
ILC_meta<- ENV[1:47,]
NZC_meta<- ENV[48:94,]

#remove columns with only zeros
ILC_notnull<-ILC[,colSums(ILC) !=0]
NZC_notnull<-NZC[,colSums(NZC) !=0]

#remove rows with only zeros, remove redundant prefixes from row.names
ILC_notnull2<-ILC_notnull[rowSums(ILC_notnull) !=0,]
NZC_notnull2<-NZC_notnull[rowSums(NZC_notnull) !=0,]

##############################################################
#Calculate some stats for original matrices
##############################################################

#get total number of taxa per site
ILC_S<-specnumber(ILC_notnull2)
NZC_S<-specnumber(NZC_notnull2)
#get largest library size with most reads
ILC_readmax<-max(rowSums(ILC_notnull2))
NZC_readmax<-max(rowSums(NZC_notnull2))
#get smallest library size with least reads
ILC_readmin<-min(rowSums(ILC_notnull2))
NZC_readmin<-min(rowSums(NZC_notnull2))
#get largest library richness with most taxa
ILC_taxmax<-max(ILC_S)
NZC_taxmax<-max(NZC_S)
#calculate 15th percentile for rrarefy function
ILC_15percentile<-quantile(rowSums(ILC_notnull2), prob=0.15)
NZC_15percentile<-quantile(rowSums(NZC_notnull2), prob=0.15)

###################################################################
##### Plot rarefaction curves
###################################################################

#print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile
# pdf("ILC_rarecurve.pdf")
# ILC_rarecurveout<-rarecurve_ESV(ILC_notnull2, sample=ILC_15percentile, step=100, col="blue4", cex=0.6)
# dev.off()
# pdf("NZC_rarecurve.pdf")
# NZC_rarecurveout<-rarecurve_ESV(NZC_notnull2, sample=NZC_15percentile, step=100, col="red4", cex=0.6)
# dev.off()

###################################################################
##### Rarefy the dataset down to the 15th percentile
###################################################################

#Rarefy original matrix down to 15th percentile library size to normalize read depth across samples
ILC_df<-rrarefy(ILC_notnull2,sample=ILC_15percentile)
NZC_df<-rrarefy(NZC_notnull2,sample=NZC_15percentile)

###################################################################
##### Calculate some stats
###################################################################

#Calculate some stats for rarefied matrices
#Get total number of species per site
ILC_S2<-specnumber(ILC_df)
NZC_S2<-specnumber(NZC_df)
#Get largest library size with most reads
ILC_readmax2<-max(rowSums(ILC_df))
NZC_readmax2<-max(rowSums(NZC_df))
#Get smallest library size with least reads
ILC_readmin2<-min(rowSums(ILC_df))
NZC_readmin2<-min(rowSums(NZC_df))
#Get largest library richness with most taxa
ILC_taxmax2<-max(ILC_S2)
NZC_taxmax2<-max(NZC_S2)

###################################################################
##### Print stome stats before and after rarefaction
###################################################################

#Compare original and rarefied stats
cat("Smallest library size (reads)")
cat("\nOriginal ILC: ", ILC_readmin)
cat("\nRarefied ILC: ", ILC_readmin2)
cat("\nOriginal NZC: ", NZC_readmin)
cat("\nRarefied NZC: ", NZC_readmin2)
cat("\nLargest library richness (ESVs)")
cat("\nOriginal ILC: ", ILC_taxmax)
cat("\nRarefied ILC: ", ILC_taxmax2)
cat("\nOriginal NZC: ", NZC_taxmax)
cat("\nRarefied NZC: ", NZC_taxmax2)

###################################################################
##### Transform to presence absence matrix
###################################################################

#Transform into a presence-absence matrix
ILC_df[ILC_df>0] <-1
NZC_df[NZC_df>0] <-1

###################################################################
##### Alpha diversity richness boxplots
###################################################################

library(ggplot2)
library(reshape)

#turn matrix into dataframe so it can be manipulated
ILC_df_df<-data.frame(ILC_df)
dimshape<-dim(ILC_df_df)
#create new empty df to summarize data from ILC_df
ILC_table<-data.frame()[1:dimshape[[1]],]
row.names(ILC_table)<-row.names(ILC_df_df)
dim(ILC_table)

#add data to new df
ILC_table$ESVrichness<-rowSums(ILC_df_df)
ILC_table$expt<-c(rep("1E",each=23),
                 rep("3E",each=24))
ILC_table$layer<-c(rep("B",each=8),rep("M",each=8),rep("O",each=7),
                   rep("B",each=8),rep("M",each=8),rep("O",each=8))
ILC_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55")

#combine expt and layer columns
ILC_table$expt_layer<-paste(ILC_table$expt, ILC_table$layer, sep="_")

#just keep the richness, gridcoord, expt_layer columns
ILC_table<-ILC_table[c("gridcoord","expt","expt_layer","ESVrichness")]
head(ILC_table)
dim(ILC_table)

#indicate factor ordering, show pooled across layers (no power to distinguish among layers)
ILC_table$expt_layer<-factor(ILC_table$expt_layer,
                             levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                              ordered=TRUE)
box<-ggplot(ILC_table, aes(x=expt,y=ESVrichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (ESVs)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))

pdf("ILC_alpha_boxplot.pdf")
box
dev.off()

#turn matrix into dataframe so it can be manipulated
NZC_df_df<-data.frame(NZC_df)
dimshape2<-dim(NZC_df_df)
#create new empty df to summarize data from ILC_df
NZC_table<-data.frame()[1:dimshape2[[1]],]
row.names(NZC_table)<-row.names(NZC_df_df)
row.names(NZC_table)
dim(NZC_table)

#add data to new df
NZC_table$ESVrichness<-rowSums(NZC_df_df)
NZC_table$expt<-c(rep("1E",each=24),
                  rep("3E",each=23))
NZC_table$layer<-c(rep("B",each=8),rep("M",each=8),rep("O",each=8),
                   rep("B",each=8),rep("M",each=8),rep("O",each=7))
NZC_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","43","51","55")

#combine expt and layer columns
NZC_table$expt_layer<-paste(NZC_table$expt, NZC_table$layer, sep="_")

#just keep the richness, gridcoord, expt_layer columns
NZC_table<-NZC_table[c("gridcoord","expt","expt_layer","ESVrichness")]
head(NZC_table)
dim(NZC_table)

#indicate factor ordering, show pooled across layers (no power to distinguish among layers)
NZC_table$expt_layer<-factor(NZC_table$expt_layer,
                             levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                             ordered=TRUE)
box2<-ggplot(NZC_table, aes(x=expt,y=ESVrichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (ESVs)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))

pdf("NZC_alpha_boxplot.pdf")
box2
dev.off()
###################################################################
##### do 3 paired-t tests 
###################################################################

#reshape ILC_table
ILC_table2<-cast(ILC_table,gridcoord ~ expt_layer, sum)
head(ILC_table2)
class(ILC_table2)
#create new df
ILC_table3<-data.frame()[1:24,]
dim(ILC_table3)

#pool data across layers
E1_B<-as.vector(ILC_table2[[2]])
E1_O<-as.vector(ILC_table2[[4]])
E1_M<-as.vector(ILC_table2[[6]])
E1<-c(E1_B,E1_O,E1_M)
ILC_table3$E1<-E1
E3_B<-as.vector(ILC_table2[[3]])
E3_O<-as.vector(ILC_table2[[5]])
E3_M<-as.vector(ILC_table2[[7]])
E3<-c(E3_B,E3_O,E3_M)
ILC_table3$E3<-E3
print(ILC_table3)
ILC_table4<-ILC_table3[-10,]
print(ILC_table4)

# turn dataframe into matrix so each column is numeric vector
DF<-data.matrix(ILC_table4, rownames.force=NA)
print(DF)

#Compare E1 and E3 (pooled across layers) ##### GOES WITH BOXPLOT #####
t.test(DF[,1],DF[,2],paired=TRUE)
#t=0.091619, df=22, p=0.9278

#reshape NZC_table
NZC_table2<-cast(NZC_table,gridcoord ~ expt_layer, sum)
head(NZC_table2)
class(NZC_table2)
#create new df
NZC_table3<-data.frame()[1:24,]
dim(NZC_table3)

#pool data across layers
N_E1_B<-as.vector(NZC_table2[[2]])
N_E1_O<-as.vector(NZC_table2[[4]])
N_E1_M<-as.vector(NZC_table2[[6]])
N_E1<-c(N_E1_B,N_E1_O,N_E1_M)
NZC_table3$N_E1<-N_E1
N_E3_B<-as.vector(NZC_table2[[3]])
N_E3_O<-as.vector(NZC_table2[[5]])
N_E3_M<-as.vector(NZC_table2[[7]])
N_E3<-c(N_E3_B,N_E3_O,N_E3_M)
NZC_table3$N_E3<-N_E3
print(NZC_table3)
NZC_table4<-NZC_table3[-13,]
print(NZC_table4)

# turn dataframe into matrix so each column is numeric vector
N_DF<-data.matrix(NZC_table4, rownames.force=NA)
print(N_DF)

#Compare E1 and E3 (pooled across layers) ##### GOES WITH BOXPLOT #####
t.test(N_DF[,1],N_DF[,2],paired=TRUE)
#t=-2.6714, df=22, p=0.0.01394

###################################################################
##### Power analysis, t test
###################################################################

library(pwr)

#Cohen 1988 t=0.2 small, t=0.5 medium, t=0.8 large effect sizes

#What sample size would we need to detect a large effect size sensu Cohen 1988
#pwr.t.test(d=0.8, sig.level=0.05, power=0.8, type="paired")
#n=14, ok so pool 1E and 3E samples n=24 (23 paired)

#If we pool data across layers, what power do we have to detect large effect if there?
#throw out unpaired val
pwr.t.test(n=23, d=0.8, sig.level=0.05, type="paired")
#power=0.963 DO T-TEST THIS WAY

#Check power for small effect, data pooled across layers #ILC
pwr.t.test(n=23, d=0.091619, sig.level=0.05, type="paired")
#power=0.07 not enough to detect small effect

#Check power for small effect, data pooled across layers #NZC
pwr.t.test(n=23, d=-2.6714, sig.level=0.05, type="paired")
#power=1 enough to detect a large effect

#Check power to detect med effect, data pooled across layers
pwr.t.test(n=23, d=0.5, sig.level=0.05, type="paired")
#power=0.63, eh
##### BASICALLY WE ONLY HAVE ENOUGH POWER TO DETECT LARGE EFFECTS #####

#Would we have power to detect med effect if we pooled data across two sites?
pwr.t.test(n=48, d=0.5, sig.level=0.05, type="paired")
#power=0.92, yes

#Would we have power to detect small effect if we pooled data across two sites?
pwr.t.test(n=48, d=0.2, sig.level=0.05, type="paired")
#power=0.27, no

###################################################################
##### Plot sample size curves for detecting significant effects w/ t test
###################################################################

# range of t
r <- seq(.1,.9,.1)
nr <- length(r)

# power values
p <- seq(.4,.9,.1)
np <- length(p)

# obtain sample sizes
samsize <- array(numeric(nr*np), dim=c(nr,np))
for (i in 1:np){
  for (j in 1:nr){
    result <- pwr.t.test(n = NULL, d= r[j], sig.level = .05, power = p[i])
    samsize[j,i] <- ceiling(result$n)}}

# set up graph
xrange <- range(r)
yrange <- round(range(samsize))

#get ggplot colors
ggplotColors<-function(g){
  d<-360/g
  h<-cumsum(c(15,rep(d,g-1)))
  hcl(h=h,c=100,l=65)
}
colors<-ggplotColors(6)

pdf("PowerAnalysisPairedTtest.pdf")
plot(xrange, yrange, type="n", xlab="Effect size (t)", ylab="Sample Size (n)" ,
     cex.lab=1.2,cex.axis=1.2)

# add power curves
for (i in 1:np){
  lines(r, samsize[,i], type="l", lwd=2, col=colors[i])
}

# add annotation (grid lines, title, legend)
abline(h=23, lty=2, col="black")
title("Power analysis, paired T-test, Sig=0.05")
legend("topright", title="Power", as.character(p), fill=colors, bty="n")
dev.off()

###################################################################
##### Create distance matrix for beta diversity analyses
###################################################################

ILC_sor<-vegdist(ILC_df, "bray", binary=TRUE)
NZC_sor<-vegdist(NZC_df, "bray", binary=TRUE)

###################################################################
##### Richness bar plot
###################################################################

#read infile prepared by python script
I_richness<-read.csv(file="ILC_richness.csv", head=TRUE, row.names=1)
I_richness

N_richness<-read.csv(file="NZC_richness.csv", head=TRUE, row.names=1)
N_richness

library(ggplot2)
library(reshape)

#reshape ILC data for ggplot
I_richness2<-melt(I_richness, 
                id.vars=c("rank"))
I_richness2

Ig <- ggplot(I_richness2, aes(fill=variable,y=value,x=rank)) +
  geom_bar(position="dodge",stat="identity") +
  theme_bw() +
  ggtitle("ILC45") +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        plot.title=element_text(size=20, hjust=0.5)) +
  xlab("Rank") +
  ylab("Richness (log10)") +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_discrete(name="Extractions",
                      breaks=c("X1E","X3E"),
                      labels=c("1","3"))
  
pdf("ILC_richness.pdf")
Ig
dev.off()

#reshape NZC data for ggplot
N_richness2<-melt(N_richness, 
                  id.vars=c("rank"))
N_richness2

Ng <- ggplot(N_richness2, aes(fill=variable,y=value,x=rank)) +
  geom_bar(position="dodge",stat="identity") +
  theme_bw() +
  ggtitle("NZC85") +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        plot.title=element_text(size=20, hjust=0.5)) +
  xlab("Rank") +
  ylab("Richness (log10)") +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_discrete(name="Extractions",
                      breaks=c("X1E","X3E"),
                      labels=c("1","3"))

pdf("NZC_richness.pdf")
Ng
dev.off()
###################################################################
##### UPGMA clustering
#################################################s##################

#create distance matrix based on P-A data using jaccard dissimilarity
# pdf("ILC_Sorensen_upgma.pdf")
# #print(row.names(ILC_df))
# print(row.names(ILC_df)<-gsub("ILC45_","",row.names(ILC_df)))
# ILC_jac<-vegdist(ILC_df, "bray",binary=TRUE) #already done above
# dend<-as.dendrogram(hclust(ILC_jac, method="average")) %>%
#  set("by_labels_branches_lwd",
#    c("1C3E_O_11","1C3E_O_15","1C3E_O_23","1C3E_O_31","1C3E_O_35","1C3E_O_43","1C3E_O_51","1C3E_O_55",
#    "1C3E_B_11","1C3E_B_15","1C3E_B_23","1C3E_B_31","1C3E_B_35","1C3E_B_43","1C3E_B_51","1C3E_B_55",
#    "1C3E_M_11","1C3E_M_15","1C3E_M_23","1C3E_M_31","1C3E_M_35","1C3E_M_43","1C3E_M_51","1C3E_M_55"),
#    TF_values=c(3,Inf)) %>%
#  set("by_labels_branches_col",
#      c("1C3E_B_11","1C3E_B_15","1C3E_B_23","1C3E_B_31","1C3E_B_35","1C3E_B_43","1C3E_B_51","1C3E_B_55",
#         "1C1E_B_11","1C1E_B_15","1C1E_B_23","1C1E_B_31","1C1E_B_35","1C1E_B_43","1C1E_B_51","1C1E_B_55"),
#       TF_values=c("blue4",Inf)) %>%
#   set("by_labels_branches_col",
#       c("1C3E_O_11","1C3E_O_15","1C3E_O_23","1C3E_O_31","1C3E_O_35","1C3E_O_43","1C3E_O_51","1C3E_O_55",
#         "1C1E_O_11","1C1E_O_23","1C1E_O_31","1C1E_O_35","1C1E_O_43","1C1E_O_51","1C1E_O_55"),
#       TF_values=c("green4",Inf)) %>%
#   set("by_labels_branches_col",
#       c("1C3E_M_11","1C3E_M_15","1C3E_M_23","1C3E_M_31","1C3E_M_35","1C3E_M_43","1C3E_M_51","1C3E_M_55",
#         "1C1E_M_11","1C1E_M_15","1C1E_M_23","1C1E_M_31","1C1E_M_35","1C1E_M_43","1C1E_M_51","1C1E_M_55"),
#       TF_values=c("red4",Inf))
# #labels(dend)
# labels(dend)<-c(3,5,7,7,1,1,6,6,3,5,4,4,
#                 8,2,2,8,8,4,4,1,1,7,7,3,
#                 3,2,2,5,5,6,6,8,1,1,4,4,
#                 3,3,2,5,5,8,6,6,8,7,7)
# #labels(dend)
# par(mar=c(5,3,4,2), cex=0.8)
# plot(dend,
#      horiz=T,
#      main="ILC",
#      xlab="Sorensen dissimilarity",
# #     ylab="ILC",
#      sub=NA,
#      axes=F,
#      cex.main=2,
#      cex.lab=2,
# #     cex.axis=2,
#      pt.cex=1,
#      cex=1)
# axis(side=1, at=seq(0,1.0,0.25), cex.axis=1.5)
# dev.off()
#ILC_jac
#class(ILC_jac)

#############################################################
#Assess dispersion (variance) using PCoA
#############################################################

BD_Ext<-betadisper(ILC_sor, factor(ILC_meta$Extractions))

pdf("ILC_betadisp.pdf")
par(mfrow=c(2,2))
#Principal Coordinates Analysis used to visualize dissimilarities = Multidimensional Scaling, MDS
#shows that the distances of the group cetnroids same for 1E and 3E, 
#see 3 clusters that correspond to 3 layers
plot(BD_Ext, main="ILC beta dispersion, Sorensen", sub="", col=colors)
op<-par(cex=0.6)
legend("topright", c("1 DNA extraction","3 DNA extractions"), text.col=colors, bty="n")
#shows distance to group centroid, same for 1E and 3E
boxplot(BD_Ext)
dev.off()

#Test the null hypothesis that the gruop mean distances are the same, need a high F to reject

#method 1: ANOVA
#anova(BD_Ext)

#method2: permutation test
#permutest(BD_Ext, pairwise=TRUE, permutations=99)

#method3: Tukey's Honest Significant Difference method
#TukeyHSD(BD_Ext)

#significant difference found among Layers B/O/M FOR COMPARISON ONLY
#BD_Layer<-betadisper(ILC_jac, factor(ILC_meta$Layer))
#head(BD_Layer)
#anova(BD_Layer)
#permutest(BD_Layer, pairwise=TRUE, permutations=99)
#TukeyHSD(BD_Layer)
#plot(BD_Layer)
#boxplot(BD_Layer)

#pdf("BetaDispersion_PCoA.pdf")
#par(mfrow=c(2,2))
#plot(BD_Ext)
#plot(BD_Layer)
#dev.off()

#pdf("BetaDispersion_boxplot.pdf")
#par(mfrow=c(2,2))
#boxplot(BD_Ext)
#boxplot(BD_Layer)
#dev.off()

N_BD_Ext<-betadisper(NZC_sor, factor(NZC_meta$Extractions))

pdf("NZC_betadisp.pdf")
par(mfrow=c(2,2))
#Principal Coordinates Analysis used to visualize dissimilarities = Multidimensional Scaling, MDS
#shows that the distances of the group cetnroids same for 1E and 3E, 
#see 3 clusters that correspond to 3 layers
plot(N_BD_Ext,main="NZC beta dispersion, Sorensen",sub="",col=colors)
op<-par(cex=0.6)
legend("topright", c("1 DNA extraction","3 DNA extractions"), text.col=colors, bty="n")
#shows distance to group centroid, same for 1E and 3E
boxplot(N_BD_Ext)
dev.off()

#Test the null hypothesis that the gruop mean distances are the same, need a high F to reject

#method 1: ANOVA
#anova(BD_Ext)

#method2: permutation test
#permutest(BD_Ext, pairwise=TRUE, permutations=99)

#method3: Tukey's Honest Significant Difference method
#TukeyHSD(BD_Ext)

#significant difference found among Layers B/O/M FOR COMPARISON ONLY
#BD_Layer<-betadisper(ILC_jac, factor(ILC_meta$Layer))
#head(BD_Layer)
#anova(BD_Layer)
#permutest(BD_Layer, pairwise=TRUE, permutations=99)
#TukeyHSD(BD_Layer)
#plot(BD_Layer)
#boxplot(BD_Layer)

#pdf("BetaDispersion_PCoA.pdf")
#par(mfrow=c(2,2))
#plot(BD_Ext)
#plot(BD_Layer)
#dev.off()

#pdf("BetaDispersion_boxplot.pdf")
#par(mfrow=c(2,2))
#boxplot(BD_Ext)
#boxplot(BD_Layer)
#dev.off()
###################################################################
##### Pearson correlation coefficient of ESVs P-A
###################################################################

pdf("ILC_pearson_ESVs.pdf")

#compute a correlation matrix, use P-A matrix, transposed to look for correlations among samples
#clean up sample names
print(row.names(ILC_df)<-gsub("ILC45_", "", row.names(ILC_df)))
print(row.names(ILC_df)<-gsub("1C3E_","",row.names(ILC_df)))
print(row.names(ILC_df)<-gsub("1C1E_","",row.names(ILC_df)))
ILC_df_t<-t(ILC_df)
print(ILC_df_t)[1:5,1:5]
print(colnames(ILC_df_t))

#create matrix of just 1E and one for 3E
ILC_df_t_1E<-ILC_df_t[,1:23]
print(colnames(ILC_df_t_1E))
ILC_df_t_3E<-ILC_df_t[,24:47]
print(colnames(ILC_df_t_3E))

ILC_cor<-cor(ILC_df_t_1E,ILC_df_t_3E,method="pearson",use="pairwise.complete.obs")
ILC_cor
x.scale <- list(cex=0.5, alternating=1, col='black') 
y.scale <- list(cex=0.5, alternating=1, col='black') 
bp1 <- levelplot(ILC_cor,xlab="1 DNA extraction",ylab="3 DNA extractions",cex=0.08,    
                 at=do.breaks(c(-1.01,1.01),101),
                 scales=list(x=list(rot=90)),
                 colorkey=list(space="top"), 
                 col.regions=colorRampPalette(c("blue4","white","green4")), 
                 main=list(label=paste(substitute("Island Lake ESVs"), "\n", 
                                       substitute("1 vs 3 DNA extractions"), sep=""),cex=1))
print(bp1)
dev.off()

pdf("NZC_pearson.pdf")

#compute a correlation matrix, use P-A matrix, transposed to look for correlations among samples
#clean up sample names
print(row.names(NZC_df))
print(row.names(NZC_df)<-gsub("NZC85_", "", row.names(NZC_df)))
print(row.names(NZC_df)<-gsub("1C3E_","",row.names(NZC_df)))
print(row.names(NZC_df)<-gsub("1C1E_","",row.names(NZC_df)))
NZC_df_t<-t(NZC_df)
print(NZC_df_t)[1:5,1:5]
print(colnames(NZC_df_t))

#create matrix of just 1E and one for 3E
NZC_df_t_1E<-NZC_df_t[,1:24]
print(colnames(NZC_df_t_1E))
NZC_df_t_3E<-NZC_df_t[,25:47]
print(colnames(NZC_df_t_3E))

NZC_cor<-cor(NZC_df_t_1E,NZC_df_t_3E,method="pearson",use="pairwise.complete.obs")
NZC_cor
x.scale <- list(cex=0.5, alternating=1, col='black')
y.scale <- list(cex=0.5, alternating=1, col='black')
bp2 <- levelplot(NZC_cor,xlab="1 DNA extraction",ylab="3 DNA extractions",cex=0.08,
                at=do.breaks(c(-1.01,1.01),101),
                scales=list(x=list(rot=90)),
                colorkey=list(space="top"),
                col.regions=colorRampPalette(c("blue4","white","green4")),
                main=list(label=paste(substitute("Nimitz ESVs"), "\n",
                                      substitute("1 vs 3 DNA extractions"), sep=""),cex=1))
print(bp2)
dev.off()

##########################################################
##### Indicator species analysis
##########################################################

library(indicspecies)
library(reshape2)

#turn data frame into table
class(ILC_df)
head(ILC_df)[,1:5]
ILC_df2<-data.frame(ILC_df)
#rownames(benthic)<-benthic[,1]
class(ILC_df2)
head(ILC_df2)[,1:5]
row.names(ILC_df2)

#group1 = 1E, group2 = 3E
groups<-c(rep(1,23), 
          rep(2,24))
indval=multipatt(ILC_df2, groups, control=how(nperm=999))
summary(indval)

#turn data frame into table
class(NZC_df)
head(NZC_df)[,1:5]
NZC_df2<-data.frame(NZC_df)
#rownames(benthic)<-benthic[,1]
class(NZC_df2)
head(NZC_df2)[,1:5]
row.names(NZC_df2)

#group1 = 1E, group2 = 3E
groups<-c(rep(1,24), 
          rep(2,23))
indval=multipatt(NZC_df2, groups, control=how(nperm=999))
summary(indval)
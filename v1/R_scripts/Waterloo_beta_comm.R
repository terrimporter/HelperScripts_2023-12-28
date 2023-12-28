library(vegan)
require(graphics)
require(dendextend)
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
A<-read.csv(file="matrix.csv", head=TRUE, row.names=1)

#remove columns with only zeros
A_notnull<-A[,colSums(A) !=0]

#remove rows with only zeros, remove redundant prefixes from row.names
A_notnull2<-A_notnull[rowSums(A_notnull) !=0,]

#Calculate some stats for original matrices
#get total number of taxa per site
A_S<- specnumber(A_notnull2)

#get largest library size with most reads
A_readmax<-max(rowSums(A_notnull2))

#get smallest library size with least reads
A_readmin<-min(rowSums(A_notnull2))

#get largest library richness with most taxa
A_taxmax<-max(A_S)

#calculate 15th percentile for rrarefy function
A_15percentile<-quantile(rowSums(A_notnull2), prob=0.15)
#print(rowSums(A_notnull2))
#print(A_15percentile)

###################################################################
##### Plot rarefaction curves
###################################################################

#print rarefaction curves for each sample to assess read coverage per sample, use 15th percentile
#pdf("rarecurve.pdf")
#A_rarecurveout<-rarecurve_ESV(A_notnull2, sample=A_15percentile, step=100, col="blue4", cex=0.6)
#dev.off()


###################################################################

#Rarefy original matrix down to smallest library size to normalize read depth across samples to 15th percentile
A_df<-rrarefy(A_notnull2,sample=A_15percentile)

#Calculate some stats for rarefied matrices
#Get total number of species per site
A_S2<-specnumber(A_df)

#Get largest library size with most reads
A_readmax2<-max(rowSums(A_df))

#Get smallest library size with least reads
A_readmin2<-min(rowSums(A_df))

#Get largest library richness with most taxa
A_taxmax2<-max(A_S2)

###################################################################
##### Print stome stats before and aftr rarefaction
###################################################################

#Compare original and rarefied stats
cat("Smallest library size (reads)")
cat("\nOriginal: ", A_readmin)
cat("\nRarefied: ", A_readmin2)
cat("\nLargest library richness (ESVs)")
cat("\nOriginal: ", A_taxmax)
cat("\nRarefied: ", A_taxmax2)
cat("\n15th percentile library size for rrarefy:")
cat(A_15percentile)

###################################################################
##### Do scree plots
###################################################################

#Scree plots to determine number of dimensions to use for NMDS

#Transform into a presence-absence matrix
A_df[A_df>0] <-1

#pdf("scree.pdf")
#nmds.scree(A_notnull2)
#dev.off()

###################################################################
##### Do NMDS
###################################################################

#Transform into a presence-absence matrix, alreaady done before doing scree
#A_df[A_df>0] <-1

#Do 2 dimensional NMDS, no environment file for now 
pdf("nmds2.pdf")
A_nmds2<-metaMDS(A_df, k=3,trymax=100)

#Try to fit environmental variables
ef<-envfit(A_nmds2, ENV, permu=999)
ef

A_nmds2
shrink<-FALSE
 
#create list of colors for each point
print(row.names(A_df))
par(mfrow=c(2,2))

#1-14 urban, 16-24 natural
cols<-c(rep("black",each=20), 
        rep("red4", each=20),
        rep("black",each=6),
        rep("red4",each=18),
        rep("black",each=35))

#grab sites & scores
A_spp.sc <- scores(A_nmds2, display = "species", shrink = shrink)
A_site.sc <- scores(A_nmds2, display = "sites") 
 
# work out the ranges NMDS1 vs NMDS2
A_ylim12 <- range(A_spp.sc[,2], A_site.sc[,2], na.rm=TRUE)
A_xlim12 <- range(A_spp.sc[,1], A_site.sc[,1], na.rm=TRUE)
 
## set-up the plotting region and add points NMDS1 vs NMDS2
plot(A_site.sc, xlim = A_xlim12, ylim = A_ylim12, type = "n", 
      xlab = "NMDS1", ylab = "NMDS2", main="Waterloo Benthic Samples", cex.axis=1.5, cex.lab=1.5, cex.main=2) 

#Nothing to plot here, fit is very bad
#plot(df, p.max=0.1)
text(A_site.sc[,1:2], row.names(A_df),cex=1,col=cols)
text(-1.5, 1, "Stress = 0.165", cex=1)
text(-1.5, 0.9, expression("Linear"~R^2~"= 0.874"), cex=1)

# work out the ranges NMDS2 vs NMDS3
A_ylim23 <- range(A_spp.sc[,3], A_site.sc[,3], na.rm=TRUE)
A_xlim23 <- range(A_spp.sc[,2], A_site.sc[,2], na.rm=TRUE)

## set-up the plotting region and add points
plot(A_site.sc, xlim = A_xlim23, ylim = A_ylim23, type = "n", 
     xlab = "NMDS2", ylab = "NMDS3", main="Waterloo Benthic Samples", cex.axis=1.5, cex.lab=1.5, cex.main=2) 
text(A_site.sc[,2:3], row.names(A_df),cex=1,col=cols)
text(-1.5, 1, "Stress = 0.165", cex=1)
text(-1.5, 0.9, expression("Linear"~R^2~"= 0.874"), cex=1)

# work out the ranges NMDS3 vs NMDS1
A_ylim31 <- range(A_spp.sc[,1], A_site.sc[,1], na.rm=TRUE)
A_xlim31 <- range(A_spp.sc[,3], A_site.sc[,3], na.rm=TRUE)

## set-up the plotting region and add points
plot(A_site.sc, xlim = A_xlim31, ylim = A_ylim31, type = "n", 
     xlab = "NMDS3", ylab = "NMDS1", main="Waterloo Benthic Samples", cex.axis=1.5, cex.lab=1.5, cex.main=2) 
text(A_site.sc[,3], A_site.sc[,1], row.names(A_df),cex=1,col=cols)
text(-1.5, 1, "Stress = 0.166", cex=1)
text(-1.5, 0.9, expression("Linear"~R^2~"= 0.874"), cex=1)

#legend
#legend("topright", 
#        legend=c("1 DNA extraction", "3 DNA extractions", "Bryophyte", "Organic", "Mineral"), 
#        col=c("black","black","blue4","green4","red4"), 
#        pch=c(1,16,16,16,16), 
#        cex=1, pt.cex=2,
#        bty="o")
# 
dev.off()

###################################################################
##### Shephards curve and goodness of fit calcs
###################################################################

#Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
#pdf("stressplot.pdf")
#stressplot(A_nmds2)
#gof <-goodness(A_nmds2)
#gof
#plot(A_nmds2, display = "sites", type="n")
#points(A_nmds2, display="sites",cex=2*gof/mean(gof))
#dev.off()

###################################################################
##### Do beta diversity variance tests
###################################################################

#Read in metadata
ENV<-read.table(file="Metadata.csv", sep=",",head=TRUE, row.names=1)

#create distance matrix based on P-A data using jaccard dissimilarity
pdf("jaccard_upgma.pdf")
A_jac<-vegdist(A_df, "bray",binary=TRUE)
dend<-as.dendrogram(hclust(A_jac, method="average"))
par(mar=c(5,3,4,2), cex=0.8)
plot(dend,
     horiz=T,
     main="Waterloo Benthic Samples",
     xlab="Sorensen dissimilarity",
     sub=NA,
     axes=F,
     cex.main=2,
     cex.lab=2,
     cex=1)
axis(side=1, at=seq(0,1.0,0.25), cex.axis=1.5)
dev.off()

#Assess dispersion (variance) using ANOVA #### NEED METADATA #### check for differences among Settings (urban vs natural)
BD_Setting<-betadisper(A_jac, factor(ENV$Setting))

#Principal Coordinates Analysis used to visualize dissimilarities = Multidimensional Scaling, MDS
#shows that the distances of the group cetnroids
#plot(BD_Setting)

#shows distance to group centroid, same for 1E and 3E
#boxplot(BD_Setting)

#Test the null hypothesis that the gruop mean distances are the same, need a high F to reject
#method 1: ANOVA
anova(BD_Setting)

#method2: permutation test
permutest(BD_Setting, pairwise=TRUE, permutations=99)

#method3: Tukey's Honest Significant Difference method
TukeyHSD(BD_Setting)

#check for significant differences among sampling month (May vs August)
BD_Month<-betadisper(A_jac, factor(ENV$Month))
head(BD_Month)
anova(BD_Month)
permutest(BD_Month, pairwise=TRUE, permutations=99)
TukeyHSD(BD_Month)
#plot(BD_Month)
#boxplot(BD_Month)

#check for significant differences among sampling sites (1-riffle downstream, 2-pond, 3- riffle upstream)
BD_Rep<-betadisper(A_jac, factor(ENV$Rep))
head(BD_Rep)
anova(BD_Rep)
permutest(BD_Rep, pairwise=TRUE, permutations=99)
TukeyHSD(BD_Rep)
#plot(BD_Rep)
#boxplot(BD_Rep)

pdf("BetaDispersion_PCoA.pdf")
par(mfrow=c(2,2))
plot(BD_Setting)
plot(BD_Month)
plot(BD_Rep)
dev.off()

pdf("BetaDispersion_boxplot.pdf")
par(mfrow=c(2,2))
boxplot(BD_Setting)
boxplot(BD_Month)
boxplot(BD_Rep)
dev.off()

###################################################################
##### Pearson correlation coefficient of ESVs P-A
###################################################################

#pdf("pearson_ESVs.pdf")

#compute a correlation matrix, use P-A matrix, transposed to look for correlations among samples
#clean up sample names
#A_df_t<-t(A_df)

#A_cor<-cor(A_df_t,A_df_t,method="pearson",use="pairwise.complete.obs")

#x.scale <- list(cex=0.5, alternating=1, col='black') 
#y.scale <- list(cex=0.5, alternating=1, col='black') 
#bp1 <- levelplot(A_cor,xlab="",ylab="",cex=0.08,    
#                 at=do.breaks(c(-1.01,1.01),101),
#                 scales=list(x=list(rot=90)),
#                 colorkey=list(space="top"), 
#                 col.regions=colorRampPalette(c("blue4","white","green4")), 
#                 main=list(label=paste(substitute("Waterloo Benthic Samples"), "\n", 
#                                       substitute(""), sep=""),cex=1))
#print(bp1)
#dev.off()

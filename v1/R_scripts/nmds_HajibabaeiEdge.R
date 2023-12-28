library(vegan)
require(graphics)
require(dendextend)
library(ggplot2)
library(ggpubr)
library(corrplot)

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
##### EDIT VEGAN RARECURVE FUNCTIONS TO ADJUST XLAB AND YLAB #######
###################################################################

rarecurve_taxa<-function (x, step, sample,main="",xlab="Reads", ylab="Taxa",label=TRUE, col, lty, ...){
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
A<-read.csv(file="matrix2.csv", head=TRUE, row.names=1)

#remove last column of sums
A_df<-A[,-1]

#read infile metadata
ENV_A<-read.table(file="metadata2.csv", sep=",", head=TRUE, row.names=1)

#remove columns with only zeros
A_notnull<-A_df[,colSums(A_df) !=0]

#remove rows with only zeros, remove redundant prefixes from row.names
A_notnull2<-A_notnull[rowSums(A_notnull) !=0,]

#figure out 15th percentile sample size
A_15percentile<-quantile(rowSums(A_notnull2), prob=0.15)

###################################################################
##### Plot rarefaction curves
###################################################################

#print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile
pdf("ILC_rarecurve_taxa.pdf")
A_rarecurveout<-rarecurve_taxa(A_notnull2, sample=A_15percentile, step=100, col="blue4", cex=0.6)
dev.off()
### waaaay too many taxa excluded so work with proportions instead of rarefying

###################################################################
# convert to proportions
###################################################################

A_prop<-cbind(A_df[1], prop.table(as.matrix(A_df[-1]), margin = 1))

###################################################################
##### Create distance matrix for beta diversity analyses ESVs only
###################################################################

bray<-vegdist(A_prop, "bray", binary=FALSE)

###################################################################
##### UPGMA clustering
#################################################s##################

#create distance matrix based on P-A data using sorensen dissimilarity
pdf("Bray_upgma.pdf")
print(row.names(A_prop))
samples<-row.names(A_prop)
set1<-samples[1:26]
set2<-samples[27:43]
set3<-samples[44:52]
dend<-as.dendrogram(hclust(bray, method="average")) %>%
 set("by_labels_branches_col",
     set1,
      TF_values=c("blue4",Inf)) %>%
  set("by_labels_branches_col",
      set2,
      TF_values=c("green4",Inf)) %>%
  set("by_labels_branches_col",
      set3,
      TF_values=c("red4",Inf))
par(mar=c(5,3,4,2), cex=0.8)
plot(dend,
     horiz=T,
     main="",
     xlab="Bray-Curtis dissimilarity",
     ylab="",
     sub=NA,
     axes=F,
     cex.main=2,
     cex.lab=2,
     cex=1)
axis(side=1, at=seq(0,1.0,0.25), cex.axis=1.5)
dev.off()

#############################################################
#Assess dispersion (variance) using PCoA
#############################################################

BD_Dat<-betadisper(bray, factor(ENV_A$Date))
BD_Loc<-betadisper(bray, factor(ENV_A$Subtype))

pdf("betadisp_bray.pdf")
par(mfrow=c(2,2))
#Principal Coordinates Analysis used to visualize dissimilarities = Multidimensional Scaling, MDS
#shows that the distances of the group cetnroids same for 1E and 3E, 
#see 3 clusters that correspond to 3 layers
Dat_PCoA<-plot(BD_Dat, main="Beta dispersion, Bray-Curtis", sub="", col=c('#F8766D','#00C094','#53B400'))
legend("topright", c("6/27/13","7/4/13","7/9/13"), text.col=c('#F8766D','#00C094','#53B400'), bty="n")

Lay_PCoA<-plot(BD_Loc, main="Beta dispersion, Bray-Curtis", sub="", col=c('#00B6EB','#A58AFF','#C49A00'))
legend("topright", c("EC", "MC", "EC-Outfall"), text.col=c('#00B6EB','#A58AFF','#C49A00'), bty="n")

#shows distance to group centroid, same for 1E and 3E
Dat_box<-boxplot(BD_Dat, col=c('#F8766D','#00C094','#53B400'))
Loc_box<-boxplot(BD_Loc, col=c('#00B6EB','#A58AFF','#C49A00'),las=2)
dev.off()

#Test the null hypothesis that the gruop mean distances are the same, need a high F to reject

#method 1: ANOVA
anova(BD_Dat)
anova(BD_Loc)

#method2: permutation test
permutest(BD_Dat, pairwise=TRUE, permutations=99)
permutest(BD_Loc, pairwise=TRUE, permutations=99) ## MC & EC-outfall p.obs=0.026, permuted p=0.05

#method3: Tukey's Honest Significant Difference method
TukeyHSD(BD_Dat)
TukeyHSD(BD_Loc)

###################################################################
##### Do scree plots
###################################################################

#Scree plots to determine number of dimensions to use for NMDS
#pdf("scree.pdf")
#nmds.scree(A_prop)
#dev.off()

######################################################
# NMDS
######################################################

#Do 2 dimensional NMDS, no environment file for now 
pdf("nmds2_bray.pdf")
nmds2<-metaMDS(A_prop, k=3,trymax=100)
shrink<-FALSE
#nmds when k=2 = 0.19 *not very satisfactory
#nmds when k=3 = 0.137 *better

#create list of colors for each point
#1E = red, 3E = black
print(row.names(A_prop))

#Color indicates date
#6/27/13 red4, 7/4/13 green4, 7/9/13 blue4
cols<-c(rep("red4",each=7), 
        rep("green4", each=10), 
        rep("blue4", each=9),
        rep("red4",each=7), 
        rep("green4", each=1),
        rep("blue4", each=9),
        rep("red4",each=4),
        rep("green4", each=3),
        rep("blue4", each=2))
#Symbol Location Subtype
#open circles (1) EC, open squares (0) MC, filled circless (16) EC-outflow
pch<-c(rep(1,each=26),rep(0,each=17), rep(16,each=9))

#grab sites & scores
spp.sc <- scores(nmds2, display = "species", shrink = shrink)
site.sc <- scores(nmds2, display = "sites")

## work out the ranges
ylim <- range(spp.sc[,2], site.sc[,2], na.rm=TRUE)
xlim <- range(spp.sc[,1], site.sc[,1], na.rm=TRUE)

## set-up the plotting region and add points
plot(site.sc, xlim = xlim, ylim = ylim, type = "n",
     xlab = "NMDS1", ylab = "NMDS2", main="", cex.axis=1.5, cex.lab=1.5, cex.main=2)
points(site.sc, col=cols, pch=pch, cex=4)
#text(-1.75,-1,"Stress=0.19")

#legend
legend("topright",
       legend=c("6/27/13 drier", "7/4/13 drier", "7/9/13 wetter", "EC", "MC","EC-outflow"),
       text.col=c("red4","green4","blue4","black","black","black"),
       pch=c(NA,NA,NA,1,0,16),
       cex=1, pt.cex=2,
       bty="o")
dev.off()

###################################################################
##### Shephards curve and goodness of fit calcs
###################################################################

#Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

library(vegan)
require(graphics)
require(dendextend)

###################################################################
##### EDIT VEGAN RARECURVE FUNCTIONS TO ADJUST XLAB AND YLAB #######
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

rarecurve_genus<-function (x, step, sample,main="",xlab="Reads", ylab="Genera",label=TRUE, col, lty, ...){
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
rarecurve_family<-function (x, step, sample,main="",xlab="Reads", ylab="Families",label=TRUE, col, lty, ...){
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
B<-read.csv(file="1C1E_1C3E_genus_matrix.csv", head=TRUE, row.names=1)
C<-read.csv(file="1C1E_1C3E_family_matrix.csv", head=TRUE, row.names=1)

#split out ILC and NZC datasets
ILC_A<-A[1:47,]
NZC_A<-A[48:94,]
ILC_B<-B[1:43,]
NZC_B<-B[44:87,]
ILC_C<-C[1:44,]
NZC_C<-C[45:91,]

#read infile metadata
ENV_A<-read.table(file="file_ESVs.env", sep=",", head=TRUE, row.names=1)
ENV_B<-read.table(file="file_genera.env", sep=",", head=TRUE, row.names=1)
ENV_C<-read.table(file="file_families.env", sep=",", head=TRUE, row.names=1)

#create metadata file for ILC and NZC
ILC_A_meta<- ENV_A[1:47,]
NZC_A_meta<- ENV_A[48:94,]
ILC_B_meta<-ENV_B[1:43,]
NZC_B_meta<-ENV_B[44:87,]
ILC_C_meta<-ENV_C[1:44,]
NZC_C_meta<-ENV_C[45:91,]

#remove columns with only zeros
ILC_A_notnull<-ILC_A[,colSums(ILC_A) !=0]
NZC_A_notnull<-NZC_A[,colSums(NZC_A) !=0]
ILC_B_notnull<-ILC_B[,colSums(ILC_B) !=0]
NZC_B_notnull<-NZC_B[,colSums(NZC_B) !=0]
ILC_C_notnull<-ILC_C[,colSums(ILC_C) !=0]
NZC_C_notnull<-NZC_C[,colSums(NZC_C) !=0]

#remove rows with only zeros, remove redundant prefixes from row.names
ILC_A_notnull2<-ILC_A_notnull[rowSums(ILC_A_notnull) !=0,]
NZC_A_notnull2<-NZC_A_notnull[rowSums(NZC_A_notnull) !=0,]
ILC_B_notnull2<-ILC_B_notnull[rowSums(ILC_B_notnull) !=0,]
NZC_B_notnull2<-NZC_B_notnull[rowSums(NZC_B_notnull) !=0,]
ILC_C_notnull2<-ILC_C_notnull[rowSums(ILC_C_notnull) !=0,]
NZC_C_notnull2<-NZC_C_notnull[rowSums(NZC_C_notnull) !=0,]

##############################################################
#Calculate some stats for original matrices
##############################################################

#get total number of taxa per site
ILC_A_S<-specnumber(ILC_A_notnull2)
NZC_A_S<-specnumber(NZC_A_notnull2)
ILC_B_S<-specnumber(ILC_B_notnull2)
NZC_B_S<-specnumber(NZC_B_notnull2)
ILC_C_S<-specnumber(ILC_C_notnull2)
NZC_C_S<-specnumber(NZC_C_notnull2)

#get largest library size with most reads
ILC_A_readmax<-max(rowSums(ILC_A_notnull2))
NZC_A_readmax<-max(rowSums(NZC_A_notnull2))
ILC_B_readmax<-max(rowSums(ILC_B_notnull2))
NZC_B_readmax<-max(rowSums(NZC_B_notnull2))
ILC_C_readmax<-max(rowSums(ILC_C_notnull2))
NZC_C_readmax<-max(rowSums(NZC_C_notnull2))

#get smallest library size with least reads
ILC_A_readmin<-min(rowSums(ILC_A_notnull2))
NZC_A_readmin<-min(rowSums(NZC_A_notnull2))
ILC_B_readmin<-min(rowSums(ILC_B_notnull2))
NZC_B_readmin<-min(rowSums(NZC_B_notnull2))
ILC_C_readmin<-min(rowSums(ILC_C_notnull2))
NZC_C_readmin<-min(rowSums(NZC_C_notnull2))

#get largest library richness with most taxa
ILC_A_taxmax<-max(ILC_A_S)
NZC_A_taxmax<-max(NZC_A_S)
ILC_B_taxmax<-max(ILC_B_S)
NZC_B_taxmax<-max(NZC_B_S)
ILC_C_taxmax<-max(ILC_C_S)
NZC_C_taxmax<-max(NZC_C_S)

#calculate 15th percentile for rrarefy function
ILC_A_15percentile<-quantile(rowSums(ILC_A_notnull2), prob=0.15)
NZC_A_15percentile<-quantile(rowSums(NZC_A_notnull2), prob=0.15)
ILC_B_15percentile<-quantile(rowSums(ILC_B_notnull2), prob=0.15)
NZC_B_15percentile<-quantile(rowSums(NZC_B_notnull2), prob=0.15)
ILC_C_15percentile<-quantile(rowSums(ILC_C_notnull2), prob=0.15)
NZC_C_15percentile<-quantile(rowSums(NZC_C_notnull2), prob=0.15)

###################################################################
##### Plot rarefaction curves
###################################################################

#print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile
pdf("ILC_rarecurve_ESVs.pdf")
ILC_A_rarecurveout<-rarecurve_ESV(ILC_A_notnull2, sample=ILC_A_15percentile, step=100, col="blue4", cex=0.6)
dev.off()
pdf("NZC_rarecurve_ESVs.pdf")
NZC_A_rarecurveout<-rarecurve_ESV(NZC_A_notnull2, sample=NZC_A_15percentile, step=100, col="red4", cex=0.6)
dev.off()
pdf("ILC_rarecurve_genera.pdf")
ILC_B_rarecurveout<-rarecurve_genus(ILC_B_notnull2, sample=ILC_B_15percentile, step=100, col="blue4", cex=0.6)
dev.off()
pdf("NZC_rarecurve_genera.pdf")
NZC_B_rarecurveout<-rarecurve_genus(NZC_B_notnull2, sample=NZC_B_15percentile, step=100, col="red4", cex=0.6)
dev.off()
pdf("ILC_rarecurve_families.pdf")
ILC_C_rarecurveout<-rarecurve_family(ILC_C_notnull2, sample=ILC_C_15percentile, step=100, col="blue4", cex=0.6)
dev.off()
pdf("NZC_rarecurve_families.pdf")
NZC_C_rarecurveout<-rarecurve_family(NZC_C_notnull2, sample=NZC_C_15percentile, step=100, col="red4", cex=0.6)
dev.off()

###################################################################
##### Rarefy the dataset down to the 15th percentile
###################################################################

#Rarefy original matrix down to 15th percentile library size to normalize read depth across samples
ILC_A_df<-rrarefy(ILC_A_notnull2,sample=ILC_A_15percentile)
NZC_A_df<-rrarefy(NZC_A_notnull2,sample=NZC_A_15percentile)

### Only do this for ESVs becaue too much missing data for genus and family ranks ###
# ILC_B_df<-rrarefy(ILC_B_notnull2,sample=ILC_B_15percentile)
# NZC_B_df<-rrarefy(NZC_B_notnull2,sample=NZC_B_15percentile)
# ILC_C_df<-rrarefy(ILC_C_notnull2,sample=ILC_C_15percentile)
# NZC_C_df<-rrarefy(NZC_C_notnull2,sample=NZC_C_15percentile)

###################################################################
##### Calculate some stats
###################################################################

#Calculate some stats for rarefied matrices
#Get total number of species per site
ILC_A_S2<-specnumber(ILC_A_df)
NZC_A_S2<-specnumber(NZC_A_df)
ILC_B_S2<-specnumber(ILC_B_df)
NZC_B_S2<-specnumber(NZC_B_df)
ILC_C_S2<-specnumber(ILC_C_df)
NZC_C_S2<-specnumber(NZC_C_df)

#Get largest library size with most reads
ILC_A_readmax2<-max(rowSums(ILC_A_df))
NZC_A_readmax2<-max(rowSums(NZC_A_df))
ILC_B_readmax2<-max(rowSums(ILC_B_df))
NZC_B_readmax2<-max(rowSums(NZC_B_df))
ILC_C_readmax2<-max(rowSums(ILC_C_df))
NZC_C_readmax2<-max(rowSums(NZC_C_df))

#Get smallest library size with least reads
ILC_A_readmin2<-min(rowSums(ILC_A_df))
NZC_A_readmin2<-min(rowSums(NZC_A_df))
ILC_B_readmin2<-min(rowSums(ILC_B_df))
NZC_B_readmin2<-min(rowSums(NZC_B_df))
ILC_C_readmin2<-min(rowSums(ILC_C_df))
NZC_C_readmin2<-min(rowSums(NZC_C_df))

#Get largest library richness with most taxa
ILC_A_taxmax2<-max(ILC_A_S2)
NZC_A_taxmax2<-max(NZC_A_S2)
ILC_B_taxmax2<-max(ILC_B_S2)
NZC_B_taxmax2<-max(NZC_B_S2)
ILC_C_taxmax2<-max(ILC_C_S2)
NZC_C_taxmax2<-max(NZC_C_S2)

###################################################################
##### Print stome stats before and after rarefaction
###################################################################

#Compare original and rarefied stats
# cat("Smallest library size (reads)")
# cat("\nOriginal ILC: ", ILC_readmin)
# cat("\nRarefied ILC: ", ILC_readmin2)
# cat("\nOriginal NZC: ", NZC_readmin)
# cat("\nRarefied NZC: ", NZC_readmin2)
# cat("\nLargest library richness (ESVs)")
# cat("\nOriginal ILC: ", ILC_taxmax)
# cat("\nRarefied ILC: ", ILC_taxmax2)
# cat("\nOriginal NZC: ", NZC_taxmax)
# cat("\nRarefied NZC: ", NZC_taxmax2)

###################################################################
##### Transform to presence absence matrix
###################################################################

#Transform rarified matrix into a presence-absence matrix for Beta div calcns
ILC_A_df[ILC_A_df>0] <-1
NZC_A_df[NZC_A_df>0] <-1
ILC_B_df[ILC_B_df>0] <-1
NZC_B_df[NZC_B_df>0] <-1
ILC_C_df[ILC_C_df>0] <-1
NZC_C_df[NZC_C_df>0] <-1

#Transform non-rarified matrix into a presence-absence matrix for Alpha div calcns
ILC_A_notnull2[ILC_A_notnull2>0] <- 1
NZC_A_notnull2[NZC_A_notnull2>0] <- 1
ILC_B_notnull2[ILC_B_notnull2>0] <- 1
NZC_B_notnull2[NZC_B_notnull2>0] <- 1
ILC_C_notnull2[ILC_C_notnull2>0] <- 1
NZC_C_notnull2[NZC_C_notnull2>0] <- 1

###################################################################
##### Alpha diversity richness boxplots
###################################################################

library(ggplot2)
library(reshape)

#turn matrix into dataframe so it can be manipulated
ILC_A_df_df<-data.frame(ILC_A_notnull2)
dimshape_A<-dim(ILC_A_df_df)
ILC_B_df_df<-data.frame(ILC_B_notnull2)
dimshape_B<-dim(ILC_B_df_df)
ILC_C_df_df<-data.frame(ILC_C_notnull2)
dimshape_C<-dim(ILC_C_df_df)

#create new empty df to summarize data from ILC_df
ILC_A_table<-data.frame()[1:dimshape_A[[1]],]
row.names(ILC_A_table)<-row.names(ILC_A_df_df)
dim(ILC_A_table)
ILC_B_table<-data.frame()[1:dimshape_B[[1]],]
row.names(ILC_B_table)<-row.names(ILC_B_df_df)
dim(ILC_B_table)
ILC_C_table<-data.frame()[1:dimshape_C[[1]],]
row.names(ILC_C_table)<-row.names(ILC_C_df_df)
dim(ILC_C_table)

#add data to new df
ILC_A_table$ESVrichness<-rowSums(ILC_A_df_df)
ILC_A_table$expt<-c(rep("1E",each=23),
                 rep("3E",each=24))
ILC_A_table$layer<-c(rep("B",each=8),rep("M",each=8),rep("O",each=7),
                   rep("B",each=8),rep("M",each=8),rep("O",each=8))
ILC_A_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55")
ILC_B_table$GenusRichness<-rowSums(ILC_B_df_df)
ILC_B_table$expt<-c(rep("1E",each=20),
                    rep("3E",each=23))
ILC_B_table$layer<-c(rep("B",each=8),rep("M",each=5),rep("O",each=7),
                     rep("B",each=8),rep("M",each=7),rep("O",each=8))
ILC_B_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                         "15","35","43","51","55",
                         "11","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","43","51","55",
                         "11","15","23","31","35","43","51","55")
ILC_C_table$FamilyRichness<-rowSums(ILC_C_df_df)
ILC_C_table$expt<-c(rep("1E",each=21),
                    rep("3E",each=23))
ILC_C_table$layer<-c(rep("B",each=8),rep("M",each=6),rep("O",each=7),
                     rep("B",each=8),rep("M",each=7),rep("O",each=8))
ILC_C_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                         "11","15","35","43","51","55",
                         "11","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","43","51","55",
                         "11","15","23","31","35","43","51","55")

#combine expt and layer columns
ILC_A_table$expt_layer<-paste(ILC_A_table$expt, ILC_A_table$layer, sep="_")
ILC_B_table$expt_layer<-paste(ILC_B_table$expt, ILC_B_table$layer, sep="_")
ILC_C_table$expt_layer<-paste(ILC_C_table$expt, ILC_C_table$layer, sep="_")

#just keep the richness, gridcoord, expt_layer columns
ILC_A_table<-ILC_A_table[c("gridcoord","expt","expt_layer","ESVrichness")]
head(ILC_A_table)
dim(ILC_A_table)
ILC_B_table<-ILC_B_table[c("gridcoord","expt","expt_layer","GenusRichness")]
head(ILC_B_table)
dim(ILC_B_table)
ILC_C_table<-ILC_C_table[c("gridcoord","expt","expt_layer","FamilyRichness")]
head(ILC_C_table)
dim(ILC_C_table)

#indicate factor ordering, show pooled across layers (no power to distinguish among layers)
ILC_A_table$expt_layer<-factor(ILC_A_table$expt_layer,
                             levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                              ordered=TRUE)
box_A<-ggplot(ILC_A_table, aes(x=expt,y=ESVrichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (ESVs)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))
pdf("ILC_alpha_boxplot_ESVs.pdf")
box_A
dev.off()
ILC_B_table$expt_layer<-factor(ILC_B_table$expt_layer,
                               levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                               ordered=TRUE)
box_B<-ggplot(ILC_B_table, aes(x=expt,y=GenusRichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (Genera)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))

pdf("ILC_alpha_boxplot_genera.pdf")
box_B
dev.off()
ILC_C_table$expt_layer<-factor(ILC_C_table$expt_layer,
                               levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                               ordered=TRUE)
box_C<-ggplot(ILC_C_table, aes(x=expt,y=FamilyRichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (Families)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))
pdf("ILC_alpha_boxplot_families.pdf")
box_C
dev.off()

#turn matrix into dataframe so it can be manipulated
NZC_A_df_df<-data.frame(NZC_A_notnull2)
dimshape2_A<-dim(NZC_A_df_df)
NZC_B_df_df<-data.frame(NZC_B_notnull2)
dimshape2_B<-dim(NZC_B_df_df)
NZC_C_df_df<-data.frame(NZC_C_notnull2)
dimshape2_C<-dim(NZC_C_df_df)

#create new empty df to summarize data from NZC_df
NZC_A_table<-data.frame()[1:dimshape2_A[[1]],]
row.names(NZC_A_table)<-row.names(NZC_A_df_df)
row.names(NZC_A_table)
dim(NZC_A_table)
NZC_B_table<-data.frame()[1:dimshape2_B[[1]],]
row.names(NZC_B_table)<-row.names(NZC_B_df_df)
row.names(NZC_B_table)
dim(NZC_B_table)
NZC_C_table<-data.frame()[1:dimshape2_C[[1]],]
row.names(NZC_C_table)<-row.names(NZC_C_df_df)
row.names(NZC_C_table)
dim(NZC_C_table)

#add data to new df
NZC_A_table$ESVrichness<-rowSums(NZC_A_df_df)
NZC_A_table$expt<-c(rep("1E",each=24),
                  rep("3E",each=23))
NZC_A_table$layer<-c(rep("B",each=8),rep("M",each=8),rep("O",each=8),
                   rep("B",each=8),rep("M",each=8),rep("O",each=7))
NZC_A_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","35","43","51","55",
                       "11","15","23","31","43","51","55")
NZC_B_table$GenusRichness<-rowSums(NZC_B_df_df)
print(NZC_B_table$GenusRichness)
NZC_B_table$expt<-c(rep("1E",each=22),
                    rep("3E",each=22))
print(NZC_B_table$expt)
NZC_B_table$layer<-c(rep("B",each=8),rep("M",each=6),rep("O",each=8),
                     rep("B",each=8),rep("M",each=7),rep("O",each=7))
NZC_B_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                         "11","23","31","35","43","51",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51",
                         "11","15","23","31","43","51","55")
NZC_C_table$FamilyRichness<-rowSums(NZC_C_df_df)
NZC_C_table$expt<-c(rep("1E",each=24),
                    rep("3E",each=23))
NZC_C_table$layer<-c(rep("B",each=8),rep("M",each=8),rep("O",each=8),
                     rep("B",each=8),rep("M",each=8),rep("O",each=7))
NZC_C_table$gridcoord<-c("11","15","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","35","43","51","55",
                         "11","15","23","31","43","51","55")

#combine expt and layer columns
NZC_A_table$expt_layer<-paste(NZC_A_table$expt, NZC_A_table$layer, sep="_")
NZC_B_table$expt_layer<-paste(NZC_B_table$expt, NZC_B_table$layer, sep="_")
NZC_C_table$expt_layer<-paste(NZC_C_table$expt, NZC_C_table$layer, sep="_")

#just keep the richness, gridcoord, expt_layer columns
NZC_A_table<-NZC_A_table[c("gridcoord","expt","expt_layer","ESVrichness")]
head(NZC_A_table)
dim(NZC_A_table)
NZC_B_table<-NZC_B_table[c("gridcoord","expt","expt_layer","GenusRichness")]
head(NZC_B_table)
dim(NZC_B_table)
NZC_C_table<-NZC_C_table[c("gridcoord","expt","expt_layer","FamilyRichness")]
head(NZC_C_table)
dim(NZC_C_table)

#indicate factor ordering, show pooled across layers (no power to distinguish among layers)
NZC_A_table$expt_layer<-factor(NZC_A_table$expt_layer,
                             levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                             ordered=TRUE)
box2_A<-ggplot(NZC_A_table, aes(x=expt,y=ESVrichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (ESVs)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))
pdf("NZC_alpha_boxplot_ESVs.pdf")
box2_A
dev.off()
head(NZC_B_table)
NZC_B_table$expt_layer<-factor(NZC_B_table$expt_layer,
                               levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                               ordered=TRUE)
box2_B<-ggplot(NZC_B_table, aes(x=expt,y=GenusRichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (Genera)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))
pdf("NZC_alpha_boxplot_genera.pdf")
box2_B
dev.off()
NZC_C_table$expt_layer<-factor(NZC_C_table$expt_layer,
                               levels=c("1E_B","3E_B","1E_O","3E_O","1E_M","3E_M"),
                               ordered=TRUE)
box2_C<-ggplot(NZC_C_table, aes(x=expt,y=FamilyRichness, fill=expt)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of DNA extractions") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  ylab("Richness (Families)") +
  theme(text=element_text(size=16)) +
  theme(axis.text.y=element_text(size=14)) +
  scale_x_discrete(labels=c("1","3"))
pdf("NZC_alpha_boxplot_families.pdf")
box2_C
dev.off()

###################################################################
##### do 3 paired-t tests 
###################################################################

#reshape ILC_table
ILC_A_table2<-cast(ILC_A_table,gridcoord ~ expt_layer, sum)
head(ILC_A_table2)
class(ILC_A_table2)
ILC_B_table2<-cast(ILC_B_table,gridcoord ~ expt_layer, sum)
head(ILC_B_table2)
class(ILC_B_table2)
ILC_C_table2<-cast(ILC_C_table,gridcoord ~ expt_layer, sum)
head(ILC_C_table2)
class(ILC_C_table2)

#create new df
ILC_A_table3<-data.frame()[1:24,]
dim(ILC_A_table3)
ILC_B_table3<-data.frame()[1:24,]
dim(ILC_B_table3)
ILC_C_table3<-data.frame()[1:24,]
dim(ILC_C_table3)

#pool data across layers
E1_A_B<-as.vector(ILC_A_table2[[2]])
E1_A_O<-as.vector(ILC_A_table2[[4]])
E1_A_M<-as.vector(ILC_A_table2[[6]])
E1_A<-c(E1_A_B,E1_A_O,E1_A_M)
ILC_A_table3$E1<-E1_A
E3_A_B<-as.vector(ILC_A_table2[[3]])
E3_A_O<-as.vector(ILC_A_table2[[5]])
E3_A_M<-as.vector(ILC_A_table2[[7]])
E3_A<-c(E3_A_B,E3_A_O,E3_A_M)
ILC_A_table3$E3<-E3_A
print(ILC_A_table3)
class(ILC_A_table3)
row_sub_A=apply(ILC_A_table3, 1, function(row) all(row !=0))
ILC_A_table4<-ILC_A_table3[row_sub_A,]
#ILC_A_table4<-ILC_A_table3[-10,]
print(ILC_A_table4)
E1_B_B<-as.vector(ILC_B_table2[[2]])
E1_B_O<-as.vector(ILC_B_table2[[4]])
E1_B_M<-as.vector(ILC_B_table2[[6]])
E1_B<-c(E1_B_B,E1_B_O,E1_B_M)
ILC_B_table3$E1<-E1_B
E3_B_B<-as.vector(ILC_B_table2[[3]])
E3_B_O<-as.vector(ILC_B_table2[[5]])
E3_B_M<-as.vector(ILC_B_table2[[7]])
E3_B<-c(E3_B_B,E3_B_O,E3_B_M)
ILC_B_table3$E3<-E3_B
print(ILC_B_table3)
row_sub_B=apply(ILC_B_table3, 1, function(row) all(row !=0))
ILC_B_table4<-ILC_B_table3[row_sub_B,]
#ILC_B_table4<-ILC_B_table3[-10,]
print(ILC_B_table4)
E1_C_B<-as.vector(ILC_C_table2[[2]])
E1_C_O<-as.vector(ILC_C_table2[[4]])
E1_C_M<-as.vector(ILC_C_table2[[6]])
E1_C<-c(E1_C_B,E1_C_O,E1_C_M)
ILC_C_table3$E1<-E1_C
E3_C_B<-as.vector(ILC_C_table2[[3]])
E3_C_O<-as.vector(ILC_C_table2[[5]])
E3_C_M<-as.vector(ILC_C_table2[[7]])
E3_C<-c(E3_C_B,E3_C_O,E3_C_M)
ILC_C_table3$E3<-E3_C
print(ILC_C_table3)
row_sub_C=apply(ILC_C_table3, 1, function(row) all(row !=0))
ILC_C_table4<-ILC_C_table3[row_sub_C,]
#ILC_C_table4<-ILC_C_table3[-10,]
print(ILC_C_table4)

# turn dataframe into matrix so each column is numeric vector
DF_A<-data.matrix(ILC_A_table4, rownames.force=NA)
print(DF_A)
DF_B<-data.matrix(ILC_B_table4, rownames.force=NA)
print(DF_B)
DF_C<-data.matrix(ILC_C_table4, rownames.force=NA)
print(DF_C)

#Compare E1 and E3 (pooled across layers) ##### GOES WITH BOXPLOT #####
t.test(DF_A[,1],DF_A[,2],paired=TRUE)
#t=0.15098, df=22, p=0.8814
t.test(DF_B[,1],DF_B[,2],paired=TRUE)
#t=-0.32769, df=18, p=0.7469
t.test(DF_C[,1],DF_C[,2],paired=TRUE)
#t=-1.2811, df=19, p=0.2156 ### large effect size, NO sig diff, report this one

#reshape NZC_table
NZC_A_table2<-cast(NZC_A_table,gridcoord ~ expt_layer, sum)
head(NZC_A_table2)
class(NZC_A_table2)
NZC_B_table2<-cast(NZC_B_table,gridcoord ~ expt_layer, sum)
head(NZC_B_table2)
class(NZC_B_table2)
NZC_C_table2<-cast(NZC_C_table,gridcoord ~ expt_layer, sum)
head(NZC_C_table2)
class(NZC_C_table2)

#create new df
NZC_A_table3<-data.frame()[1:24,]
dim(NZC_A_table3)
NZC_B_table3<-data.frame()[1:24,]
dim(NZC_B_table3)
NZC_C_table3<-data.frame()[1:24,]
dim(NZC_C_table3)

#pool data across layers
N_E1_A_B<-as.vector(NZC_A_table2[[2]])
N_E1_A_O<-as.vector(NZC_A_table2[[4]])
N_E1_A_M<-as.vector(NZC_A_table2[[6]])
N_E1_A<-c(N_E1_A_B,N_E1_A_O,N_E1_A_M)
NZC_A_table3$N_E1<-N_E1_A
N_E3_A_B<-as.vector(NZC_A_table2[[3]])
N_E3_A_O<-as.vector(NZC_A_table2[[5]])
N_E3_A_M<-as.vector(NZC_A_table2[[7]])
N_E3_A<-c(N_E3_A_B,N_E3_A_O,N_E3_A_M)
NZC_A_table3$N_E3<-N_E3_A
print(NZC_A_table3)
row_sub_A<-apply(NZC_A_table3, 1, function(row) all(row !=0))
NZC_A_table4<-NZC_A_table3[row_sub_A,]
#NZC_table4<-NZC_table3[-13,]
print(NZC_A_table4)
N_E1_B_B<-as.vector(NZC_B_table2[[2]])
N_E1_B_O<-as.vector(NZC_B_table2[[4]])
N_E1_B_M<-as.vector(NZC_B_table2[[6]])
N_E1_B<-c(N_E1_B_B,N_E1_B_O,N_E1_B_M)
NZC_B_table3$N_E1<-N_E1_B
N_E3_B_B<-as.vector(NZC_B_table2[[3]])
N_E3_B_O<-as.vector(NZC_B_table2[[5]])
N_E3_B_M<-as.vector(NZC_B_table2[[7]])
N_E3_B<-c(N_E3_B_B,N_E3_B_O,N_E3_B_M)
NZC_B_table3$N_E3<-N_E3_B
print(NZC_B_table3)
row_sub_B<-apply(NZC_B_table3, 1, function(row) all(row !=0))
NZC_B_table4<-NZC_B_table3[row_sub_B,]
#NZC_table4<-NZC_table3[-13,]
print(NZC_B_table4)
N_E1_C_B<-as.vector(NZC_C_table2[[2]])
N_E1_C_O<-as.vector(NZC_C_table2[[4]])
N_E1_C_M<-as.vector(NZC_C_table2[[6]])
N_E1_C<-c(N_E1_C_B,N_E1_C_O,N_E1_C_M)
NZC_C_table3$N_E1<-N_E1_C
N_E3_C_B<-as.vector(NZC_C_table2[[3]])
N_E3_C_O<-as.vector(NZC_C_table2[[5]])
N_E3_C_M<-as.vector(NZC_C_table2[[7]])
N_E3_C<-c(N_E3_C_B,N_E3_C_O,N_E3_C_M)
NZC_C_table3$N_E3<-N_E3_C
print(NZC_C_table3)
row_sub_C<-apply(NZC_C_table3, 1, function(row) all(row !=0))
NZC_C_table4<-NZC_C_table3[row_sub_C,]
#NZC_table4<-NZC_table3[-13,]
print(NZC_C_table4)

# turn dataframe into matrix so each column is numeric vector
N_DF_A<-data.matrix(NZC_A_table4, rownames.force=NA)
print(N_DF_A)
N_DF_B<-data.matrix(NZC_B_table4, rownames.force=NA)
print(N_DF_B)
N_DF_C<-data.matrix(NZC_C_table4, rownames.force=NA)
print(N_DF_C)

#Compare E1 and E3 (pooled across layers) ##### GOES WITH BOXPLOT #####
t.test(N_DF_A[,1],N_DF_A[,2],paired=TRUE)
#t=-2.7362, df=22, p=0.01206 ### large effect size, sig diff, report this one
t.test(N_DF_B[,1],N_DF_B[,2],paired=TRUE)
#t=0.19157, df=20, p=0.85
t.test(N_DF_C[,1],N_DF_C[,2],paired=TRUE)
#t=-0.48922, df=22, p=0.6295

###################################################################
##### Power analysis, t test
###################################################################

library(pwr)

#Cohen 1988 t=0.2 small, t=0.5 medium, t=0.8 large effect sizes

#What sample size would we need to detect a large effect size sensu Cohen 1988
#pwr.t.test(d=0.8, sig.level=0.05, power=0.8, type="paired")
#n=14, ok so pool 1E and 3E samples n=24 (23 paired)

#If we pool data across layers, what power do we have to detect large effect if there?
pwr.t.test(n=24, d=0.8, sig.level=0.05, type="paired")
#power=0.963 DO T-TEST THIS WAY

#Check power to detect mod effect?
pwr.t.test(n=24, d=0.5, sig.level=0.05, type="paired")
#weak power=0.65 to detect moderate effect

#Check power for small effect
pwr.t.test(n=24, d=0.2, sig.level=0.05, type="paired")
#virtually no power=0.16  to detect a small effect

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
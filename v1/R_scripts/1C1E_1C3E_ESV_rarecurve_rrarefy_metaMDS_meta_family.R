library(vegan)

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
A<-read.csv(file="1C1E_1C3E_family_matrix.csv", head=TRUE, row.names=1)
print(row.names(A))
ENV<-read.table(file="file.env", sep=",",head=TRUE, row.names=1)

#print(ENV)
#split out ILC and NZC datasets
ILC<-A[1:46,]
NZC<-A[47:93,]

#remove columns with only zeros
ILC_notnull<-ILC[,colSums(ILC) !=0]
NZC_notnull<-NZC[,colSums(NZC) !=0]

#remove rows with only zeros, remove redundant prefixes from row.names
ILC_notnull2<-ILC_notnull[rowSums(ILC_notnull) !=0,]
#print(row.names(ILC_notnull2)<-gsub("ILC45_", "", row.names(ILC_notnull2)))
NZC_notnull2<-NZC_notnull[rowSums(NZC_notnull) !=0,]
#print(row.names(NZC_notnull2)<-gsub("NZC85_", "", row.names(NZC_notnull2)))

#Calculate some stats for original matrices
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

###################################################################
##### Plot rarefaction curves
###################################################################

#print rarefaction curves for each sample to assess read coverage per sample
#pdf("ILC_rarecurve.pdf")
#ILC_rarecurveout<-rarecurve_ESV(ILC_notnull2, sample=ILC_readmin, step=100, col="blue4", cex=0.6)
#dev.off()
#pdf("NZC_rarecurve.pdf")
#NZC_rarecurveout<-rarecurve_ESV(NZC_notnull2, sample=NZC_readmin, step=100, col="red4", cex=0.6)
#dev.off()

###################################################################

#Rarefy original matrix down to smallest library size to normalize read depth across samples
ILC_df<-rrarefy(ILC_notnull2,sample=ILC_readmin)
NZC_df<-rrarefy(NZC_notnull2,sample=NZC_readmin)

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
##### Print stome stats before and aftr rarefaction
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
##### Do scree plots
###################################################################

#Scree plots to determine number of dimensions to use for NMDS
#pdf("ILC_scree.pdf")
#nmds.scree(ILC_notnull2)
#dev.off()
#pdf("NZC_scree.pdf")
#nmds.scree(NZC_notnull2)
#dev.off()

###################################################################
##### Do NMDS
###################################################################

#Transform into a presence-absence matrix
ILC_df[ILC_df>0] <-1
NZC_df[NZC_df>0] <-1

#Do 2 dimensional NMDS, no environment file for now 
#pdf("ILC_nmds2.pdf")
#ILC_nmds2<-metaMDS(ILC_df, k=2,trymax=100)
#shrink<-FALSE

#create list of colors for each point
#1E = red, 3E = black
#print(row.names(ILC_df))
#instead of points plot grid coord, but map to single digit numbers
#ILC_gridcoord=c("11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55")
#ILC_gridcoord=c("1","2","3","4","5","6","7","8",
 #           "1","2","3","4","5","6","7","8",
  #          "1","3","4","5","6","7","8",
   #         "1","2","3","4","5","6","7","8",
    #        "1","2","3","4","5","6","7","8",
     #       "1","2","3","4","5","6","7","8")
#1=regular font for 1E, 2=bold font for 3E sampels
#ILC_fonts=c(1,1,1,1,1,1,1,1,1,1,
 #       1,1,1,1,1,1,1,1,1,1,
  #      1,1,1,2,2,2,2,2,2,2,
   #     2,2,2,2,2,2,2,2,2,2,
    #    2,2,2,2,2,2,2)
#open circles (1) 1E, filled circles (19) 3E
#ILC_pch<-c(rep(1,each=23),rep(19,each=24))
#green - bryophyte, black - organic, brown - mineral
#ILC_cols<-c(rep("green4",each=8), rep("brown4", each=8), rep("black", each=7),
 #       rep("green4",each=8), rep("brown4", each=8), rep("black", each=8))

#grab sites & scores
#ILC_spp.sc <- scores(ILC_nmds2, display = "species", shrink = shrink)
#ILC_site.sc <- scores(ILC_nmds2, display = "sites") 

## work out the ranges
#ILC_ylim <- range(ILC_spp.sc[,2], ILC_site.sc[,2], na.rm=TRUE)
#ILC_xlim <- range(ILC_spp.sc[,1], ILC_site.sc[,1], na.rm=TRUE)

## set-up the plotting region and add points
#plot(ILC_site.sc, xlim = ILC_xlim, ylim = ILC_ylim, type = "n", asp = 1,
 #    xlab = "NMDS1", ylab = "NMDS2", main="ILC 1 vs 3 DNA extractions") 
#points(ILC_site.sc, col=ILC_cols, pch=ILC_pch, cex=1.2)
#text(ILC_site.sc[,1:2], ILC_gridcoord,col=ILC_cols,cex=1, font=ILC_fonts)
#text(-0.75, 1.5, "Stress = 0.08")
#text(-0.75, 1.3, expression("Linear"~R^2~"= 0.97"))

#legend
#legend("topright",legend=c("Regular = 1 DNA extraction ", "Bold = 3 DNA extractions ", "Green = Bryophyte ", "Black = Organic ", "Red = Mineral "), 
 #      text.col=c("gray48","gray48","green4","black","brown4"), text.font=c(1,2,1,1,1), cex=0.8, pt.cex=1)

#dev.off()

#pdf("NZC_nmds2.pdf")
#NZC_nmds2<-metaMDS(NZC_df,k=2,trymax=100)
#NZC_nmds2
#shrink<-FALSE

#create list of colors for each point
#1E = red, 3E = black
#print(row.names(NZC_df))

#instead of points plot grid coord, but map to single digit numbers
#NZC_gridcoord=c("11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","15","23","31","35","43","51","55",
#            "11","23","31","35","43","51","55")
#NZC_gridcoord=c("1","2","3","4","5","6","7","8",
 #           "1","2","3","4","5","6","7","8",
  #          "1","2","3","4","5","6","7","8",
   #         "1","2","3","4","5","6","7","8",
    #        "1","2","3","4","5","6","7","8",
     #       "1","3","4","5","6","7","8")
#1=regular font for 1E, 3=bold font for 3E
#NZC_fonts=c(1,1,1,1,1,1,1,1,1,1,
 #       1,1,1,1,1,1,1,1,1,1,
  #      1,1,1,1,2,2,2,2,2,2,
   #     2,2,2,2,2,2,2,2,2,2,
    #    2,2,2,2,2,2,2)

#open circles (1) 1E, filled circles (19) 3E
#NZC_pch<-c(rep(1,each=24),rep(19,each=23))
#green - bryophyte, black - organic, brown - mineral
#NZC_cols<-c(rep("green4",each=8), rep("brown4", each=8), rep("black", each=8),
 #           rep("green4",each=8), rep("brown4", each=8), rep("black", each=7))

#grab sites & scores
#NZC_spp.sc <- scores(NZC_nmds2, display = "species", shrink = shrink)
#NZC_site.sc <- scores(NZC_nmds2, display = "sites") 

## work out the ranges
#NZC_ylim <- range(NZC_spp.sc[,2], NZC_site.sc[,2], na.rm=TRUE)
#NZC_xlim <- range(NZC_spp.sc[,1], NZC_site.sc[,1], na.rm=TRUE)

## set-up the plotting region and add points
#plot(NZC_site.sc, xlim = NZC_xlim, ylim = NZC_ylim, type = "n", asp = 1,
 #    xlab = "NMDS1", ylab = "NMDS2", main="NZC 1 vs 3 DNA extractions") 
#points(NZC_site.sc, col=NZC_cols, pch=NZC_pch, cex=1.2)
#text(NZC_site.sc[,1:2], NZC_gridcoord,col=NZC_cols,cex=1, font=NZC_fonts)
#text(-1, 1, "Stress = 0.12")
#text(-1, 0.8, expression("Linear"~R^2~"= 0.94"))

#legend
#legend("topright",legend=c("Regular = 1 DNA extraction", "Bold = 3 DNA extractions", "Green = Bryophyte", "Black = Organic", "Red = Mineral"), 
#       text.col=c("gray48","gray48","green4","black","brown4"), cex=0.8, pt.cex=1)

#dev.off()

###################################################################
##### Shephards curve and goodness of fit calcs
###################################################################

#Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
#pdf("ILC_stressplot.pdf")
#stressplot(ILC_nmds2)
#gof <-goodness(ILC_nmds2)
#gof
#plot(ILC_nmds2, display = "sites", type="n")
#points(ILC_nmds2, display="sites",cex=2*gof/mean(gof))
#dev.off()
#pdf("NZC_stressplot.pdf")
#stressplot(NZC_nmds2)
#gof <-goodness(NZC_nmds2)
#gof
#plot(NZC_nmds2, display = "sites", type="n")
#points(NZC_nmds2, display="sites",cex=2*gof/mean(gof))
#dev.off()

###################################################################
##### check integrity of metadata ENV that was imported
###################################################################

#dim(ENV)
#class(ENV)
#str(ENV[,1:3])
#names(ENV)
#rownames(ENV)[1:3]

#check if rownames in A and ENV are same, no, so change them
#head(ENV)[,1:3]
#head(A)[,1:5]
#dim(A) 
#row.names(ENV) <- row.names(A)

#sanity check, ensure rownames in original full A and ENV are identical
#identical(row.names(ENV),row.names(A)[-ncol(A)])

#create metadata file for ILC and NZC
#ILC_meta<- ENV[1:47,]
#NZC_meta<- ENV[48:94,]

#library(BiodiversityR)
#library(RAM)
#packageVersion("RAM")

#need transposed ILC_df matrix as dataframe
#class(ILC_df)
#ILC.df<-as.data.frame(ILC_df)
#class(ILC.df)
#ILC.df.t<-t(ILC.df)
#head(ILC.df.t)[1:5]
#ILC.meta.div <- OTU.diversity(data=list(ILC.df.t),meta=ILC_meta)

###################################################################
##### Do beta diversity variance tests
###################################################################

#create distance matrix based on P-A data using jaccard dissimilarity
#ILC_jac<-vegdist(ILC_df, "jac",binary=TRUE)
#ILC_jac
#class(ILC_jac)

#Assess dispersion (variance) using ANOVA
#BD_Ext<-betadisper(ILC_jac, factor(ILC_meta$Extractions))
#head(MD_Ext)

#Principal Coordinates Analysis used to visualize dissimilarities = Multidimensional Scaling, MDS
#shows that the distances of the group cetnroids same for 1E and 3E, 
#see 3 clusters that correspond to 3 layers
#plot(BD_Ext)

#shows distance to group centroid, same for 1E and 3E
#boxplot(BD_Ext)

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

pdf("ILC_pearson_family.pdf")

#compute a correlation matrix, use P-A matrix, transposed to look for correlations among samples
#clean up sample names
print(row.names(ILC_df)<-gsub("ILC45_", "", row.names(ILC_df)))
print(row.names(ILC_df)<-gsub("1C3E_","",row.names(ILC_df)))
print(row.names(ILC_df)<-gsub("1C1E_","",row.names(ILC_df)))
ILC_df_t<-t(ILC_df)
print(ILC_df_t)[1:5,1:5]
print(colnames(ILC_df_t))

#create matrix of just 1E and one for 3E
ILC_df_t_1E<-ILC_df_t[,1:22]
print(colnames(ILC_df_t_1E))
ILC_df_t_3E<-ILC_df_t[,23:44]
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
                 main=list(label=paste(substitute("Island Lake Families"), "\n", 
                                       substitute("1 vs 3 DNA extractions"), sep=""),cex=1))
print(bp1)
dev.off()

#pdf("NZC_pearson.pdf")

#compute a correlation matrix, use P-A matrix, transposed to look for correlations among samples
#clean up sample names
#print(row.names(NZC_df))
#print(row.names(NZC_df)<-gsub("NZC85_", "", row.names(NZC_df)))
#print(row.names(NZC_df)<-gsub("1C3E_","",row.names(NZC_df)))
#print(row.names(NZC_df)<-gsub("1C1E_","",row.names(NZC_df)))
#NZC_df_t<-t(NZC_df)
#print(NZC_df_t)[1:5,1:5]
#print(colnames(NZC_df_t))

#create matrix of just 1E and one for 3E
#NZC_df_t_1E<-NZC_df_t[,1:24]
#print(colnames(NZC_df_t_1E))
#NZC_df_t_3E<-NZC_df_t[,25:47]
#print(colnames(NZC_df_t_3E))

#NZC_cor<-cor(NZC_df_t_1E,NZC_df_t_3E,method="pearson",use="pairwise.complete.obs")
#NZC_cor
#x.scale <- list(cex=0.5, alternating=1, col='black') 
#y.scale <- list(cex=0.5, alternating=1, col='black') 
#bp2 <- levelplot(NZC_cor,xlab="1 DNA extraction",ylab="3 DNA extractions",cex=0.08,    
#                 at=do.breaks(c(-1.01,1.01),101),
#                 scales=list(x=list(rot=90)),
#                 colorkey=list(space="top"), 
#                 col.regions=colorRampPalette(c("blue4","white","green4")), 
#                 main=list(label=paste(substitute("Nimitz ESVs"), "\n", 
#                                       substitute("1 vs 3 DNA extractions"), sep=""),cex=1))
#print(bp2)
#dev.off()

###################################################################
##### Pearson correlation coefficient of ESVs relabund
###################################################################

pdf("ILC_pearson_relabund_family.pdf")

#convert read numbers to relative abundance/sample
print(ILC_notnull2[1:5,1:5])
ILC_df.ra<-decostand(ILC_notnull2,"total")
print(ILC_df.ra[1:5,1:5])

#compute a correlation matrix, use relabund matrix, transposed to look for correlations among samples
#clean up sample names
print(row.names(ILC_df.ra)<-gsub("ILC45_", "", row.names(ILC_df.ra)))
print(row.names(ILC_df.ra)<-gsub("1C3E_","",row.names(ILC_df.ra)))
#print(row.names(ILC_df.ra)<-gsub("1C1E_","",row.names(ILC_df.ra)))
ILC_df.ra_t<-t(ILC_df.ra)
print(ILC_df.ra_t)[1:5,1:5]
print(colnames(ILC_df.ra_t))

#create matrix of just 1E and one for 3E
ILC_df.ra_t_1E<-ILC_df.ra_t[,1:22]
print(colnames(ILC_df.ra_t_1E))
print(colnames(ILC_df.ra_t_1E)<-gsub("1C1E_","",colnames(ILC_df.ra_t_1E)))
ILC_df.ra_t_3E<-ILC_df.ra_t[,23:44]
print(colnames(ILC_df.ra_t_3E))

ILC_ra_cor<-cor(ILC_df.ra_t_1E,ILC_df.ra_t_3E,method="pearson",use="pairwise.complete.obs")
ILC_ra_cor
x.scale <- list(cex=0.5, alternating=1, col='black') 
y.scale <- list(cex=0.5, alternating=1, col='black') 
bp3 <- levelplot(ILC_ra_cor,xlab="1 DNA extraction",ylab="3 DNA extractions",cex=0.08,    
                 at=do.breaks(c(-1.01,1.01),101),
                 scales=list(x=list(rot=90)),
                 colorkey=list(space="top"), 
                 col.regions=colorRampPalette(c("blue4","white","green4")), 
                 main=list(label=paste(substitute("Island Lake Families"), "\n", 
                                       substitute("1 vs 3 DNA extractions"), sep=""),cex=1))
print(bp3)
dev.off()

#pdf("NZC_pearson_relabund.pdf")

#convert read numbers to relative abundance/sample
#print(NZC_notnull2[1:5,1:5])
#NZC_df.ra<-decostand(NZC_notnull2,"total")
#print(NZC_df.ra[1:5,1:5])

#compute a correlation matrix, use P-A matrix, transposed to look for correlations among samples
#clean up sample names
#print(row.names(NZC_df.ra))
#print(row.names(NZC_df.ra)<-gsub("NZC85_", "", row.names(NZC_df.ra)))
#print(row.names(NZC_df.ra)<-gsub("1C3E_","",row.names(NZC_df.ra)))
##print(row.names(NZC_df.ra)<-gsub("1C1E_","",row.names(NZC_df.ra)))
#NZC_df.ra_t<-t(NZC_df.ra)
#print(NZC_df.ra_t)[1:5,1:5]
#print(colnames(NZC_df.ra_t))

#create matrix of just 1E and one for 3E
#NZC_df.ra_t_1E<-NZC_df.ra_t[,1:24]
#print(colnames(NZC_df.ra_t_1E))
#print(colnames(NZC_df.ra_t_1E)<-gsub("1C1E_","",colnames(NZC_df.ra_t_1E)))
#NZC_df.ra_t_3E<-NZC_df.ra_t[,25:47]
#print(colnames(NZC_df.ra_t_3E))

#NZC_ra_cor<-cor(NZC_df.ra_t_1E,NZC_df.ra_t_3E,method="pearson",use="pairwise.complete.obs")
#NZC_ra_cor
#x.scale <- list(cex=0.5, alternating=1, col='black') 
#y.scale <- list(cex=0.5, alternating=1, col='black') 
#bp4 <- levelplot(NZC_ra_cor,xlab="1 DNA extraction",ylab="3 DNA extractions",cex=0.08,    
#                 at=do.breaks(c(-1.01,1.01),101),
#                 scales=list(x=list(rot=90)),
#                 colorkey=list(space="top"), 
#                 col.regions=colorRampPalette(c("blue4","white","green4")), 
#                 main=list(label=paste(substitute("Nimitz ESVs"), "\n", 
#                                       substitute("1 vs 3 DNA extractions"), sep=""),cex=1))
#print(bp4)
#dev.off()

###################################################################
##### Do other comparison tests
###################################################################
#also no significance found among layers
permutest(MD_Ext, pairwise=TRUE, permutations=99)

#also no significant difference found
#TukeyHSD(MD_Ext)

#significant difference found among Layers B/O/M
#MD_Layer<-betadisper(ILC_jac, factor(ILC_meta$Layer))
#head(MD_Layer)
#anova(MD_Layer)

library("asbio")
#environment(mrpp)<-environment(rda)
#no significance difference found
#MD_Ext.mrpp<-with(ILC_meta, mrpp(ILC_df, Extractions))
#MD_Ext.mrpp
#suppressWarnings(suppressMessages(suppressPackageStartupMessages(mrpp(ILC_df,ILC_meta$Extractions, permutations=999, distance="bray"))))

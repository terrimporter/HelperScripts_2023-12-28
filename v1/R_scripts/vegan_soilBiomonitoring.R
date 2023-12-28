#August 9, 2013 by Terri Porter
#Script to process multiple .csv presence-absence matrices from excel (carriage returns already removed) with 
#[R] package vegan to do scree, nmds, and stressplots
#with Biomonitoring 2.0 soil dataset (singleton's removed)

library(vegan)

########## NEW FUNCTION FOR SCREE ##########

#new function edited from http://www.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm
nmds.scree<-function(x) {
	plot(rep(1,10),replicate(10,metaMDS(x, autotransform=F,k=1)$stress/100),xlim=c(1,nrow(x)),ylim=c(0,0.5),xlab="Number dimensions",ylab="Stress",main="NMDS stress plot")
	
	for(i in 1:(nrow(x)-2)) {
		points(rep(i+1,10),replicate(10,metaMDS(x, autotransform=F,k=i+1)$stress/100))
	}
}

########## EDIT VEGAN PACKAGE STRESSPLOT FUNCTION ##########

stressplot.edit<-function(object, dis, pch, p.col = "grey", l.col = "black", lwd = 2, t.cex=0.6,
    ...) {
    require(MASS) || stop("Needs MASS package")
    if (missing(dis)) 
        dis <- metaMDSredist(object)
    if (attr(dis, "Size") != nrow(object$points)) 
        stop("Dimensions do not match in ordination and dissimilarities")
    shep <- Shepard(dis, object$points)
    stress <- sum((shep$y - shep$yf)^2)/sum(shep$y^2)
    rstress <- 1 - stress
    ralscal <- cor(shep$y, shep$yf)^2
    stress <- sqrt(stress) * 100
    if (abs(stress - object$stress) > 0.001) 
        stop("Dissimilarities and ordination do not match")
    if (missing(pch)) 
        if (length(dis) > 5000) 
            pch = "."
        else pch = 1
    plot(shep, pch = pch, col = p.col, xlab = "Observed Dissimilarity", 
        ylab = "Ordination Distance", ...)
    lines(shep$x, shep$yf, type = "S", col = l.col, lwd = lwd, 
        ...)
    lab <- paste("Non-metric fit, R2 =", format(rstress, digits = 3), 
        "\nLinear fit, R2 =", format(ralscal, digits = 3))
    text(min(shep$x), 0.95 * max(shep$y), lab, pos = 4, cex=t.cex)
    invisible(shep)
}

########## OPEN FILES ##########

rbcLF<-read.csv(file="rbcLF_incidence.csv", head=TRUE, row.names=1)
rbcLR<-read.csv(file="rbcLR_incidence.csv", head=TRUE, row.names=1)
ITS1R<-read.csv(file="ITS1R_incidence.csv", head=TRUE, row.names=1)
ITS4F<-read.csv(file="ITS4F_incidence.csv", head=TRUE, row.names=1)
SSUv3<-read.csv(file="16V3_incidence.csv", head=TRUE, row.names=1)
SSUv6<-read.csv(file="16V6_incidence.csv", head=TRUE, row.names=1)

########## DO SCREE ##########

#pdf("NMDS_scree_plots.pdf")
#par(mfrow=c(3,2))

#nmds.scree(rbcLF)
#nmds.scree(rbcLR)
#nmds.scree(ITS1R)
#nmds.scree(ITS4F)
#nmds.scree(SSUv3)
#nmds.scree(SSUv6)

#dev.off()

########## DO NMDS ##########

pdf("NMDS_plots.pdf")
par(mfcol=c(3,3)) #use 3,3 instead of 3,2 to force square plots

rbcLF.m<-metaMDS(rbcLF, distance="bray", k=3, trymax=1000, maxiter=1000)
rbcLF.m
plot(rbcLF.m, main="rbcLF", type="n")
points(rbcLF.m, display= "spec", cex = 0.35, pch=1, col="grey")
text(rbcLF.m, display= "sites", cex = 0.8, col="black")

ITS1R.m<-metaMDS(ITS1R, distance="bray", k=3, trymax=1000, maxiter=1000)
ITS1R.m
plot(ITS1R.m, main="ITS1R", type="n")
points(ITS1R.m, display= "spec", cex = 0.35, pch=1, col="grey")
text(ITS1R.m, display= "sites", cex = 0.8, col="black")

SSUv3.m<-metaMDS(SSUv3, distance="bray", k=3, trymax=1000, maxiter=1000)
SSUv3.m
plot(SSUv3.m, main="16V3", type="n")
points(SSUv3.m, display= "spec", cex = 0.35, pch=1, col="grey")
text(SSUv3.m, display= "sites", cex = 0.8, col="black")

rbcLR.m<-metaMDS(rbcLR, distance="bray", k=3, trymax=1000, maxiter=1000)
rbcLR.m
plot(rbcLR.m, main="rbcLR", type="n")
points(rbcLR.m, display= "spec", cex = 0.35, pch=1, col="grey")
text(rbcLR.m, display= "sites", cex = 0.8, col="black")

ITS4F.m<-metaMDS(ITS4F, distance="bray", k=3, trymax=1000, maxiter=1000)
ITS4F.m
plot(ITS4F.m, main="ITS4F", type="n")
points(ITS4F.m, display= "spec", cex = 0.35, pch=1, col="grey")
text(ITS4F.m, display= "sites", cex = 0.8, col="black")

SSUv6.m<-metaMDS(SSUv6, distance="bray", k=3, trymax=1000, maxiter=1000)
SSUv6.m
plot(SSUv6.m, main="16V6", type="n")
points(SSUv6.m, display= "spec", cex = 0.35, pch=1, col="grey")
text(SSUv6.m, display= "sites", cex = 0.8, col="black")

dev.off()


########## DO STRESSPLOTS ##########

pdf("NMDS_stressplots.pdf")
par(mfcol=c(3,3))

stressplot.edit(rbcLF.m, main="rbcLF")
stressplot.edit(ITS1R.m, main="ITS1R")
stressplot.edit(SSUv3.m, main="16V3")
stressplot.edit(rbcLR.m, main="rbcLR")
stressplot.edit(ITS4F.m, main="ITS4F")
stressplot.edit(SSUv6.m, main="16V6")

dev.off()

########## DO GOF PLOTS ##########

pdf("NMDS_gof.pdf")
par(mfcol=c(3,3))

gof<-goodness(rbcLF.m)
gof
plot(rbcLF.m,display="sites", type="n", main="rbcLF")
points(rbcLF.m, display="sites",cex=2*gof/mean(gof))

rm(gof)

gof<-goodness(ITS1R.m)
gof
plot(ITS1R.m,display="sites", type="n", main="ITS1R")
points(ITS1R.m, display="sites",cex=2*gof/mean(gof))

rm(gof)

gof<-goodness(SSUv3.m)
gof
plot(SSUv3.m,display="sites", type="n", main="16V3")
points(SSUv3.m, display="sites",cex=2*gof/mean(gof))

rm(gof)

gof<-goodness(rbcLR.m)
gof
plot(rbcLR.m,display="sites", type="n", main="rbcLR")
points(rbcLR.m, display="sites",cex=2*gof/mean(gof))

rm(gof)

gof<-goodness(ITS4F.m)
gof
plot(ITS4F.m,display="sites", type="n", main="ITS4F")
points(ITS4F.m, display="sites",cex=2*gof/mean(gof))

rm(gof)

gof<-goodness(SSUv6.m)
gof
plot(SSUv6.m,display="sites", type="n", main="16V6")
points(SSUv6.m, display="sites",cex=2*gof/mean(gof))

dev.off()

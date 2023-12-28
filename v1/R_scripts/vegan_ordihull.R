rbcLF<-read.csv(file="incidence.csv", head=TRUE, row.names=1)
rbcLF

rbcLF.env<-read.table(file="rbcLF.env",head=TRUE)
rbcLF.env

library(vegan)

#dat <- vegdist(rbcLF, method="bray", binary=TRUE)
#dat

#pdf("0_pooled_nmds12.pdf")
mds <- metaMDS(rbcLF, distance="bray", k=3, trymax=1000, maxiter=1000) 
#plot(mds, choices=c(1,2), main="No pooling", cex.main=0.9, type="n", axes=FALSE, cex.lab=0.8)
#axis(1, cex.axis=0.65)
#axis(2, cex.axis=0.65)
#points(mds, choices=c(1,2), display="spec", cex=0.35, pch=1, col="grey")
#points(mds, choices=c(1,2), display="sites", cex=0.35, pch=16, col="black")
#ordihull(mds, choices=c(1,2), group=rbcLF.env$Site, show="PAD", label=TRUE)
#ordihull(mds, choices=c(1,2), group=rbcLF.env$Site, show="WC", label=TRUE)
#ordispider(mds, group=rbcLF.env$Site, lty=3, col="blue")
#dev.off()

#pdf("0_pooled_nmds23.pdf")
#plot(mds, choices=c(2,3), main="No pooling", cex.main=0.9, type="n", axes=FALSE, cex.lab=0.8)
#axis(1, cex.axis=0.65)
#axis(2, cex.axis=0.65)
#points(mds, choices=c(2,3), display="spec", cex=0.35, pch=1, col="grey")
#points(mds, choices=c(2,3), display="sites", cex=0.35, pch=16, col="black")
#ordihull(mds, choices=c(2,3),  group=rbcLF.env$Site, show="PAD", label=TRUE)
#ordihull(mds, choices=c(2,3), group=rbcLF.env$Site, show="WC", label=TRUE)
#dev.off()

#pdf("0_pooled_nmds13.pdf")
#plot(mds, choices=c(1,3), main="No pooling", cex.main=0.9, type="n", axes=FALSE, cex.lab=0.8)
#axis(1, cex.axis=0.65)
#axis(2, cex.axis=0.65)
#points(mds, choices=c(1,3), display="spec", cex=0.35, pch=1, col="grey")
#points(mds, choices=c(1,3), display="sites", cex=0.35, pch=16, col="black")
#ordihull(mds, choices=c(1,3), group=rbcLF.env$Site, show="PAD", label=TRUE)
#ordihull(mds, choices=c(1,3), group=rbcLF.env$Site, show="WC", label=TRUE)
#dev.off()

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

#do stressplot for nmds1v2 only
pdf("stressplot.pdf")
stressplot.edit(mds, main="rbcLF")
dev.off()


#adonis(rbcLF ~ Site, data=rbcLF.env, strata=rbcLF.env$Site, perm=1e3)
#adonis(rbcLF ~ Site, data=rbcLF.env, perm=1e3)

#adonis(rbcLF ~ Site*SubsiteRep, data=rbcLF.env, perm=1e3)

#adonis(rbcLF ~ Site*SubsiteRep*SquareRep, data=rbcLF.env, perm=1e3)

#adonis(rbcLF ~ Site*SubsiteRep*SquareRep*CoreRep, data=rbcLF.env, perm=1e3)
	
#adonis(rbcLF ~ SubsiteRep, strata=rbcLF.env$Site, data=rbcLF.env, perm=1e3)
#adonis(rbcLF ~ SquareRep, strata=rbcLF.env$Site, data=rbcLF.env, perm=1e3)

#pdf("cluster.pdf")

#par(mfrow=c(2,2))

#dis <-vegdist(rbcLF)
#clus <- hclust(dis, "single")
#plot(clus, cex=0.5)

#clus2 <- hclust(dis,"complete")
#plot(clus2, cex=0.5)

#clus3 <-hclust(dis, "average")
#plot(clus3, cex=0.3)

#dev.off()


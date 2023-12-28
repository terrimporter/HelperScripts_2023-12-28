rbcLF<-read.csv(file="rbcLF_incidence.csv", head=TRUE, row.names=1)
rbcLF

rbcLF.env<-read.table(file="rbcLF.env",head=TRUE)
rbcLF.env

library(vegan)

#dat <- vegdist(rbcLF, method="bray", binary=TRUE)
#dat

pdf("9_pooled.pdf")

mds <- metaMDS(rbcLF)
plot(mds)

ordihull(mds, group=rbcLF.env$Site, show="PAD")
ordihull(mds, group=rbcLF.env$Site, show="WC")

#ordispider(mds, group=rbcLF.env$Site, lty=3, col="blue")

dev.off()

#adonis(rbcLF ~ Site, data=rbcLF.env, strata=rbcLF.env$Site, perm=1e3)
adonis(rbcLF ~ Site, data=rbcLF.env, perm=1e3)

adonis(rbcLF ~ Site*SubsiteRep, data=rbcLF.env, perm=1e3)

#adonis(rbcLF ~ SubsiteRep, strata=rbcLF.env$Site, data=rbcLF.env, perm=1e3)

pdf("cluster.pdf")

#par(mfrow=c(2,2))

dis <-vegdist(rbcLF)
#clus <- hclust(dis, "single")
#plot(clus)

#clus2 <- hclust(dis,"complete")
#plot(clus2)

clus3 <-hclust(dis, "average")
plot(clus3, ylab="Bray-Curtis dissimilarity", main="9 pooled cores per sample")

dev.off()

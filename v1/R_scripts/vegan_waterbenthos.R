ITS4F<-read.csv(file="ITS4F.csv", head=TRUE, row.names=1)
ITS4F

library(vegan)
ITS4F.m<-metaMDS(ITS4F, distance="bray", k=2, trymax=1000, maxiter=1000)
ITS4F.m

pdf(file="ITS4F.pdf")
plot(ITS4F.m, main="ITS4F")
dev.off()

#####

ITS1R<-read.csv(file="ITS1R.csv", head=TRUE, row.names=1)
ITS1R

ITS1R.m<-metaMDS(ITS1R, distance="bray", k=2, trymax=1000, maxiter=1000)
ITS1R.m

pdf(file="ITS1R.pdf")
plot(ITS1R.m, main="ITS1R")
dev.off()

#####

rbcLF<-read.csv(file="rbcLF.csv", head=TRUE, row.names=1)
rbcLF

rbcLF.m<-metaMDS(rbcLF, distance="bray", k=2, trymax=1000, maxiter=1000)
rbcLF.m

pdf(file="rbcLF.pdf")
plot(rbcLF.m, main="rbcLF")
dev.off()

#####

rbcLR<-read.csv(file="rbcLR.csv", head=TRUE, row.names=1)
rbcLR

rbcLR.m<-metaMDS(rbcLR, distance="bray", k=2, trymax=1000, maxiter=1000)
rbcLR.m

pdf(file="rbcLR.pdf")
plot(rbcLR.m, main="rbcLR")
dev.off()

#####

SSUv6<-read.csv(file="16v6.csv", head=TRUE, row.names=1)
SSUv6

SSUv6.m<-metaMDS(SSUv6, distance="bray", k=2, trymax=1000, maxiter=1000)
SSUv6.m

pdf(file="16Sv6.pdf")
plot(SSUv6.m, main="16Sv6")
dev.off()

#####

SSUv3<-read.csv(file="16v3.csv", head=TRUE, row.names=1)
SSUv3

SSUv3.m<-metaMDS(SSUv3, distance="bray", k=2, trymax=1000, maxiter=1000)
SSUv3.m

pdf(file="16Sv3.pdf")
plot(SSUv3.m, main="16Sv3")
dev.off()

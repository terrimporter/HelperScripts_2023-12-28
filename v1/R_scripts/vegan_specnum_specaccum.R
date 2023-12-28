rbcLF<-read.csv(file="rbcLF_incidence.csv", head=TRUE, row.names=1)
rbcLF

rbcLF.env<-read.table(file="rbcLF.env",head=TRUE)
rbcLF.env

library(vegan)

specnum<-specnumber(rbcLF)
specnum

pdf("specnum_boxplot.pdf")

boxplot(specnum ~ rbcLF.env$Site, ylab="Number of OTUs")

dev.off()

t.test(specnum ~ rbcLF.env$Site)

pdf("specaccum_samplingcurve.pdf")

par(mfrow=c(2,2))

rbcLF_PAD<-rbcLF[1:8,]
rbcLF_WC<-rbcLF[9:16,]
#remove columns with only zeros

rbcLF_PAD<-rbcLF_PAD[,colSums(rbcLF_PAD) !=0]
rbcLF_WC<-rbcLF_WC[,colSums(rbcLF_WC) !=0]

plot(specaccum(rbcLF_PAD), xlab="Number of subsites", ylab="Number of OTUs", main="PAD")
plot(specaccum(rbcLF_WC), xlab="Number of subsites", ylab="Number of OTUs", main="WC")

dev.off()

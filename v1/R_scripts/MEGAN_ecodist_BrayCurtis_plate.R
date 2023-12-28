pdf(file="ecodist_NMDS_plate.pdf")
par(mfrow=c(2,2))

#plot rbcL
a<-read.table(file="rbcL_FR_normalized_family_BC.txt")
a
labels<-rownames(a)
labels
b<-as.matrix(a)
library(Matrix)
c<-tril(b)
d<-as.matrix(c)
e<-as.dist(d)
library(ecodist)


h<-nmds(e, nits=100, maxdim=3)
i<-nmds.min(h)
i

library(RColorBrewer)
warm<-brewer.pal(9,"YlOrRd")
warm8<-tail(warm,8)
warm8r<-rev(warm8)
warm8
warm8r
warm5<-warm[c(4,6,7,8,9)]
warm5r<-rev(warm5)
warm5
warm5r
cool<-brewer.pal(9,"PuBuGn")
cool<-tail(cool,8)
coolr<-rev(cool)
cool
coolr

#pch 15 filled squares for rbcLF, pch 17 filled triangles for rbcLR
#col warm for PAD, col cool for WC
#plot(i, main="Bray-Curtis disimilarity", pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17), col=c(warm,warm,cool,cool))
#only plot dimensions 2 and 3
plot(i$X2, i$X3, main="rbcL", pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17), col=c(warm8r,warm8r,coolr,coolr), xlab="NMDS2", ylab="NMDS3")

#plot ITS
a<-read.table(file="ITS_RF_normalized_family_BC.txt")
a
labels<-rownames(a)
labels
b<-as.matrix(a)
library(Matrix)
c<-tril(b)
d<-as.matrix(c)
e<-as.dist(d)

h<-nmds(e, nits=100, maxdim=4)
i<-nmds.min(h)
i

#plot(i, main="Bray-Curtis disimilarity", pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17), col=c(warm,warm,cool,cool))
#only plot dimensions 3 and 4
plot(i$X3, i$X4, , main="ITS", pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17), col=c(warm8r,warm8r,coolr,coolr), xlab="NMDS3", ylab="NMDS4")

#plot 16S
a<-read.table(file="16v3_16v6_normalized_family_BC.txt")
a
labels<-rownames(a)
labels
b<-as.matrix(a)
library(Matrix)
c<-tril(b)
d<-as.matrix(c)
e<-as.dist(d)

h<-nmds(e, nits=100, maxdim=3)
i<-nmds.min(h)
i

#pch 16 filled circles for 16v3 pch 18 filled diamonds for 16v6
#col warm for PAD, col cool for WC
#plot(i, main="Bray-Curtis disimilarity", pch=c(16,16,16,16,16,18,18,18,18,18,18,18,18,16,16,16,16,16,16,16,16,18,18,18,18,18,18,18,18), col=c(warm5,warm8,cool,cool))
#only plot dimensions 1 and 2
plot(i$X1, i$X2, main="16S", pch=c(16,16,16,16,16,18,18,18,18,18,18,18,18,16,16,16,16,16,16,16,16,18,18,18,18,18,18,18,18), col=c(warm5r,warm8r,coolr,coolr), xlab="NMDS1", ylab="NMDS2")

#create empty plot so legend prints separately
plot(1,type="n",axes=FALSE,xlab="",ylab="")

legend("topleft", legend=c("PAD1", "PAD11", "PAD14", "PAD3", "PAD33", "PAD37", "PAD38", "PAD4", "WC1", "WC2", "WC3", "WC4", "WC5", "WC6", "WC7", "WC8", "Forward primer", "Reverse primer", "16V3 primers", "16V6 primers"), ncol=2, cex=0.75, pch=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,15,17,16,18), text.col=c(warm8r,coolr,"black","black","black","black"), col=c(warm8,cool,"black","black","black","black"), bty="n")

save.image(file="NMDS_ecodist_plate.RData")
dev.off()
quit()

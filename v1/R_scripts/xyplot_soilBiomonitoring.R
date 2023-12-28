#rbcLF
rbcLF <- read.table(file="rbcLF.csv",header=TRUE,sep=",")
rbcLF

pdf(file="rarefaction.pdf")

par(mfrow=c(3,2))
par(xpd=T, mar=par()$mar+c(0,0,0,4))

library(RColorBrewer)
warm<-brewer.pal(9,"YlOrRd")
warm<-tail(warm,8)
warm<-rev(warm)
warm
cool<-brewer.pal(9,"PuBuGn")
cool<-tail(cool,8)
cool<-rev(cool)
cool

#Plot rbcLF
plot(x=rbcLF$X, y=rbcLF$PAD1, type="l", axes=FALSE, xlab="Reads", ylab="OTUs", col=warm[1], ylim=c(0,1600), xlim=c(0,16000), main="rbcLF")
lines(x=rbcLF$X, y=rbcLF$PAD11, type="l", col=warm[2])
lines(x=rbcLF$X, y=rbcLF$PAD14, type="l", col=warm[3])
lines(x=rbcLF$X, y=rbcLF$PAD3, type="l", col=warm[4])
lines(x=rbcLF$X, y=rbcLF$PAD33, type="l", col=warm[5])
lines(x=rbcLF$X, y=rbcLF$PAD37, type="l", col=warm[6])
lines(x=rbcLF$X, y=rbcLF$PAD38, type="l", col=warm[7])
lines(x=rbcLF$X, y=rbcLF$PAD4, type="l", col=warm[8])
lines(x=rbcLF$X, y=rbcLF$WC1, type="l", col=cool[1])
lines(x=rbcLF$X, y=rbcLF$WC2, type="l", col=cool[2])
lines(x=rbcLF$X, y=rbcLF$WC3, type="l", col=cool[3])
lines(x=rbcLF$X, y=rbcLF$WC4, type="l", col=cool[4])
lines(x=rbcLF$X, y=rbcLF$WC5, type="l", col=cool[5])
lines(x=rbcLF$X, y=rbcLF$WC6, type="l", col=cool[6])
lines(x=rbcLF$X, y=rbcLF$WC7, type="l", col=cool[7])
lines(x=rbcLF$X, y=rbcLF$WC8, type="l", col=cool[8])

#x-axis
xticks<-seq(0,16000, by=2000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
axis(2,at=c(0,1600), labels=c("",""), lwd.ticks=0, pos=0)
axis(2, at=seq(0,1600,by=200), lwd=0, lwd.ticks=1, pos=0, cex.axis=0.75)

#rbcLR
rbcLR <- read.table(file="rbcLR.csv",header=TRUE,sep=",")
rbcLR

#Plot rbcLR
plot(x=rbcLR$X, y=rbcLR$PAD1, type="l", axes=FALSE, xlab="Reads", ylab="OTUs", col=warm[1], ylim=c(0,1600), xlim=c(0,16000),main="rbcLR")
lines(x=rbcLR$X, y=rbcLR$PAD11, type="l", col=warm[2])
lines(x=rbcLR$X, y=rbcLR$PAD14, type="l", col=warm[3])
lines(x=rbcLR$X, y=rbcLR$PAD3, type="l", col=warm[4])
lines(x=rbcLR$X, y=rbcLR$PAD33, type="l", col=warm[5])
lines(x=rbcLR$X, y=rbcLR$PAD37, type="l", col=warm[6])
lines(x=rbcLR$X, y=rbcLR$PAD38, type="l", col=warm[7])
lines(x=rbcLR$X, y=rbcLR$PAD4, type="l", col=warm[8])
lines(x=rbcLR$X, y=rbcLR$WC1, type="l", col=cool[1])
lines(x=rbcLR$X, y=rbcLR$WC2, type="l", col=cool[2])
lines(x=rbcLR$X, y=rbcLR$WC3, type="l", col=cool[3])
lines(x=rbcLR$X, y=rbcLR$WC4, type="l", col=cool[4])
lines(x=rbcLR$X, y=rbcLR$WC5, type="l", col=cool[5])
lines(x=rbcLR$X, y=rbcLR$WC6, type="l", col=cool[6])
lines(x=rbcLR$X, y=rbcLR$WC7, type="l", col=cool[7])
lines(x=rbcLR$X, y=rbcLR$WC8, type="l", col=cool[8])

#x-axis
xticks<-seq(0,16000, by=2000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
axis(2,at=c(0,1600), labels=c("",""), lwd.ticks=0, pos=0)
axis(2, at=seq(0,1600,by=200), lwd=0, lwd.ticks=1, pos=0, cex.axis=0.75)

#ITS1R
ITS1R <- read.table(file="ITS1R.csv",header=TRUE,sep=",")
ITS1R

#Plot ITS1R
plot(x=ITS1R$X, y=ITS1R$PAD1, type="l", axes=FALSE, xlab="Reads", ylab="OTUs", col=warm[1], ylim=c(0,1800), xlim=c(0,18000), main="ITS1R")
lines(x=ITS1R$X, y=ITS1R$PAD11, type="l", col=warm[2])
lines(x=ITS1R$X, y=ITS1R$PAD14, type="l", col=warm[3])
lines(x=ITS1R$X, y=ITS1R$PAD3, type="l", col=warm[4])
lines(x=ITS1R$X, y=ITS1R$PAD33, type="l", col=warm[5])
lines(x=ITS1R$X, y=ITS1R$PAD37, type="l", col=warm[6])
lines(x=ITS1R$X, y=ITS1R$PAD38, type="l", col=warm[7])
lines(x=ITS1R$X, y=ITS1R$PAD4, type="l", col=warm[8])
lines(x=ITS1R$X, y=ITS1R$WC1, type="l", col=cool[1])
lines(x=ITS1R$X, y=ITS1R$WC2, type="l", col=cool[2])
lines(x=ITS1R$X, y=ITS1R$WC3, type="l", col=cool[3])
lines(x=ITS1R$X, y=ITS1R$WC4, type="l", col=cool[4])
lines(x=ITS1R$X, y=ITS1R$WC5, type="l", col=cool[5])
lines(x=ITS1R$X, y=ITS1R$WC6, type="l", col=cool[6])
lines(x=ITS1R$X, y=ITS1R$WC7, type="l", col=cool[7])
lines(x=ITS1R$X, y=ITS1R$WC8, type="l", col=cool[8])

#x-axis
xticks<-seq(0,18000, by=2000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
axis(2,at=c(0,1800), labels=c("",""), lwd.ticks=0, pos=0)
axis(2, at=seq(0,1800,by=200), lwd=0, lwd.ticks=1, pos=0, cex.axis=0.75)

#ITS4F
ITS4F <- read.table(file="ITS4F.csv",header=TRUE,sep=",")
ITS4F

#Plot ITS4F
plot(x=ITS4F$X, y=ITS4F$PAD1, type="l", axes=FALSE, xlab="Reads", ylab="OTUs", col=warm[1], ylim=c(0,1800), xlim=c(0,18000), main="ITS4F")
lines(x=ITS4F$X, y=ITS4F$PAD11, type="l", col=warm[2])
lines(x=ITS4F$X, y=ITS4F$PAD14, type="l", col=warm[3])
lines(x=ITS4F$X, y=ITS4F$PAD3, type="l", col=warm[4])
lines(x=ITS4F$X, y=ITS4F$PAD33, type="l", col=warm[5])
lines(x=ITS4F$X, y=ITS4F$PAD37, type="l", col=warm[6])
lines(x=ITS4F$X, y=ITS4F$PAD38, type="l", col=warm[7])
lines(x=ITS4F$X, y=ITS4F$PAD4, type="l", col=warm[8])
lines(x=ITS4F$X, y=ITS4F$WC1, type="l", col=cool[1])
lines(x=ITS4F$X, y=ITS4F$WC2, type="l", col=cool[2])
lines(x=ITS4F$X, y=ITS4F$WC3, type="l", col=cool[3])
lines(x=ITS4F$X, y=ITS4F$WC4, type="l", col=cool[4])
lines(x=ITS4F$X, y=ITS4F$WC5, type="l", col=cool[5])
lines(x=ITS4F$X, y=ITS4F$WC6, type="l", col=cool[6])
lines(x=ITS4F$X, y=ITS4F$WC7, type="l", col=cool[7])
lines(x=ITS4F$X, y=ITS4F$WC8, type="l", col=cool[8])

#x-axis
xticks<-seq(0,18000, by=2000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
axis(2,at=c(0,1800), labels=c("",""), lwd.ticks=0, pos=0)
axis(2, at=seq(0,1800,by=200), lwd=0, lwd.ticks=1, pos=0, cex.axis=0.75)

#16V3
SSUv3 <- read.table(file="16V3.csv",header=TRUE,sep=",")
SSUv3

#Plot 16V3
plot(x=SSUv3$X, y=SSUv3$PAD14, type="l", axes=FALSE, xlab="Reads", ylab="OTUs", col=warm[3], ylim=c(0,12000), xlim=c(0,70000), main="16V3")
#lines(x=SSUv3$X, y=SSUv3$PAD11, type="l", col=warm[2])
#lines(x=SSUv3$X, y=SSUv3$PAD14, type="l", col=warm[3])
#lines(x=SSUv3$X, y=SSUv3$PAD3, type="l", col=warm[4])
lines(x=SSUv3$X, y=SSUv3$PAD33, type="l", col=warm[5])
lines(x=SSUv3$X, y=SSUv3$PAD37, type="l", col=warm[6])
lines(x=SSUv3$X, y=SSUv3$PAD38, type="l", col=warm[7])
lines(x=SSUv3$X, y=SSUv3$PAD4, type="l", col=warm[8])
lines(x=SSUv3$X, y=SSUv3$WC1, type="l", col=cool[1])
lines(x=SSUv3$X, y=SSUv3$WC2, type="l", col=cool[2])
lines(x=SSUv3$X, y=SSUv3$WC3, type="l", col=cool[3])
lines(x=SSUv3$X, y=SSUv3$WC4, type="l", col=cool[4])
lines(x=SSUv3$X, y=SSUv3$WC5, type="l", col=cool[5])
lines(x=SSUv3$X, y=SSUv3$WC6, type="l", col=cool[6])
lines(x=SSUv3$X, y=SSUv3$WC7, type="l", col=cool[7])
lines(x=SSUv3$X, y=SSUv3$WC8, type="l", col=cool[8])

#x-axis
xticks<-seq(0,70000, by=10000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
axis(2,at=c(0,12000), labels=c("",""), lwd.ticks=0, pos=0)
axis(2, at=seq(0,12000,by=2000), lwd=0, lwd.ticks=1, pos=0, cex.axis=0.75)

#16V6
SSUv6 <- read.table(file="16V6.csv",header=TRUE,sep=",")
SSUv6

#Plot 16V6
plot(x=SSUv6$X, y=SSUv6$PAD1, type="l", axes=FALSE, xlab="Reads", ylab="OTUs", col=warm[1], ylim=c(0,12000), xlim=c(0,70000), main="16V6")
lines(x=SSUv6$X, y=SSUv6$PAD11, type="l", col=warm[2])
lines(x=SSUv6$X, y=SSUv6$PAD14, type="l", col=warm[3])
lines(x=SSUv6$X, y=SSUv6$PAD3, type="l", col=warm[4])
lines(x=SSUv6$X, y=SSUv6$PAD33, type="l", col=warm[5])
lines(x=SSUv6$X, y=SSUv6$PAD37, type="l", col=warm[6])
lines(x=SSUv6$X, y=SSUv6$PAD38, type="l", col=warm[7])
lines(x=SSUv6$X, y=SSUv6$PAD4, type="l", col=warm[8])
lines(x=SSUv6$X, y=SSUv6$WC1, type="l", col=cool[1])
lines(x=SSUv6$X, y=SSUv6$WC2, type="l", col=cool[2])
lines(x=SSUv6$X, y=SSUv6$WC3, type="l", col=cool[3])
lines(x=SSUv6$X, y=SSUv6$WC4, type="l", col=cool[4])
lines(x=SSUv6$X, y=SSUv6$WC5, type="l", col=cool[5])
lines(x=SSUv6$X, y=SSUv6$WC6, type="l", col=cool[6])
lines(x=SSUv6$X, y=SSUv6$WC7, type="l", col=cool[7])
lines(x=SSUv6$X, y=SSUv6$WC8, type="l", col=cool[8])

#x-axis
xticks<-seq(0,70000, by=10000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
axis(2,at=c(0,12000), labels=c("",""), lwd.ticks=0, pos=0)
axis(2, at=seq(0,12000,by=2000), lwd=0, lwd.ticks=1, pos=0, cex.axis=0.75)

legend(75000, 12000, legend=c("PAD1", "PAD11", "PAD14", "PAD3", "PAD33", "PAD37", "PAD38", "PAD4", "WC1", "WC2", "WC3", "WC4", "WC5", "WC6", "WC7", "WC8"), cex=0.5, fill=c(warm,cool))

par(mar=c(5,4,4,2)+0.1)

dev.off()

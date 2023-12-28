#cores
#cores <- read.table(file="cores.txt",header=TRUE,sep="\t")
#cores

pdf(file="rarefaction.pdf")

par(mfrow=c(2,2))
#par(xpd=T, mar=par()$mar+c(0,0,0,4))

#library(RColorBrewer)
#warm<-brewer.pal(9,"YlOrRd")
#warm<-tail(warm,8)
#warm<-rev(warm)
#warm
#cool<-brewer.pal(9,"PuBuGn")
#cool<-tail(cool,8)
#cool<-rev(cool)
#cool

#Plot cores
#plot(x=cores$numsampled, y=cores$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,2500), ylim=c(0,1500), col="red", main="Cores (N=144)")
#lines(x=cores$numsampled, y=cores$stdevplus, type="l", col="black")
#lines(x=cores$numsampled, y=cores$stdevneg, type="l", col="black")

#x-axis
#xticks<-seq(0,2500, by=500)
#axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
#yticks<-seq(0,1500,by=500)
#axis(2,at=yticks, labels=yticks, pos=0)

#legend(250,1500, legend=c("Average", "Standard deviation"), cex=0.5, fill=c("red","black"))

#par(mar=c(5,4,4,2)+0.1)

#blocks
blocks <- read.table(file="blocks.txt",header=TRUE,sep="\t")
blocks

#Plot blocks
plot(x=blocks$numsampled, y=blocks$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,2500), ylim=c(0,700), col="red", main="Blocks (N=45)")
lines(x=blocks$numsampled, y=blocks$stdevplus, type="l", col="black")
lines(x=blocks$numsampled, y=blocks$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,2500, by=500)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,700,by=100)
axis(2,at=yticks, labels=yticks, pos=0)

legend(250,600, legend=c("Average", "Standard deviation"), cex=0.5, fill=c("red","black"))

#subsites
subsites <- read.table(file="subsites.txt",header=TRUE,sep="\t")
subsites

#Plot subsites
plot(x=subsites$numsampled, y=subsites$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,7000), ylim=c(0,1500), col="red", main="Subsites (N=16)")
lines(x=subsites$numsampled, y=subsites$stdevplus, type="l", col="black")
lines(x=subsites$numsampled, y=subsites$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,7000, by=1000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,1500,by=200)
axis(2,at=yticks, labels=yticks, pos=0)

#sites
sites <- read.table(file="sites.txt",header=TRUE,sep="\t")
sites

#Plot sites
plot(x=sites$numsampled, y=sites$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,113500), ylim=c(0,11500), col="red", main="Sites (N=2)")
lines(x=sites$numsampled, y=sites$stdevplus, type="l", col="black")
lines(x=sites$numsampled, y=sites$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,113500, by=20000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,11500,by=2000)
axis(2,at=yticks, labels=yticks, pos=0)

dev.off()

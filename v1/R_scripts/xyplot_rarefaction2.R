pdf(file="rarefaction_subsites_taxa.pdf")

par(mfrow=c(2,2))

#subsites
subsites <- read.table(file="subsites.txt",header=TRUE,sep="\t")
subsites

#Plot subsites
plot(x=subsites$numsampled, y=subsites$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,35000), ylim=c(0,10500), col="red", main="Subsites (N=16)")
lines(x=subsites$numsampled, y=subsites$stdevplus, type="l", col="black")
lines(x=subsites$numsampled, y=subsites$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,35000, by=5000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,10500,by=2500)
axis(2,at=yticks, labels=yticks, pos=0)

legend(2500,10500, legend=c("Average", "Standard deviation"), cex=0.5, fill=c("red","black"))

#plants
plants <- read.table(file="plants.txt",header=TRUE,sep="\t")
plants

#Plot plants
plot(x=plants$numsampled, y=plants$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,5000), ylim=c(0,1250), col="red", main="Subsites (Viridiplantae)")
lines(x=plants$numsampled, y=plants$stdevplus, type="l", col="black")
lines(x=plants$numsampled, y=plants$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,5000, by=1000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,1250,by=250)
axis(2,at=yticks, labels=yticks, pos=0)

#fungi
fungi <- read.table(file="fungi.txt",header=TRUE,sep="\t")
fungi

#Plot fungi
plot(x=fungi$numsampled, y=fungi$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,6000), ylim=c(0,1500), col="red", main="Subsites (Fungi)")
lines(x=fungi$numsampled, y=fungi$stdevplus, type="l", col="black")
lines(x=fungi$numsampled, y=fungi$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,6000, by=1000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,1500,by=500)
axis(2,at=yticks, labels=yticks, pos=0)

#bacteria
bacteria <- read.table(file="bacteria.txt",header=TRUE,sep="\t")
bacteria

#Plot bacteria
plot(x=bacteria$numsampled, y=bacteria$mean, type="l", axes=FALSE, xlab="Reads", ylab="Rarefied OTUs", xlim=c(0,20000), ylim=c(0,7000), col="red", main="Subsites (Bacteria)")
lines(x=bacteria$numsampled, y=bacteria$stdevplus, type="l", col="black")
lines(x=bacteria$numsampled, y=bacteria$stdevneg, type="l", col="black")

#x-axis
xticks<-seq(0,20000, by=5000)
axis(1, at=xticks, labels=xticks, pos=0)

#y-axis
yticks<-seq(0,7000,by=1000)
axis(2,at=yticks, labels=yticks, pos=0)

dev.off()

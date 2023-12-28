#September 16, 2013 by Terri Porter
#Script to create MEGAN taxonomic assignment rarefaction curves for Joel

a<-read.csv(file="malaise_90_100bs.rma-chart", header=FALSE, sep="\t")
a
b<-read.csv(file="malaise_90_200bs.rma-chart", header=FALSE, sep="\t")
b
c<-read.csv(file="malaise_95_100bs.rma-chart", header=FALSE, sep="\t")
c
d<-read.csv(file="malaise_95_200bs.rma-chart", header=FALSE, sep="\t")
d
e<-read.csv(file="malaise_98_100bs.rma-chart", header=FALSE, sep="\t")
e
f<-read.csv(file="malaise_98_200bs.rma-chart", header=FALSE, sep="\t")
f

pdf(file="rarefaction.pdf")

par(mfrow=c(3,3))

plot(a$V1, a$V2, type="l", main="90% clustering", xlab="% OTUs", ylab="Number of leaves", ylim=c(0,1600), col="red")
lines(b$V1, b$V2, type="l", col="blue")

plot(c$V1, c$V2, type="l", main="95% clustering", xlab="% OTUs", ylab="Number of leaves", ylim=c(0,1600), col="red")
lines(d$V1, d$V2, type="l", col="blue")

plot(e$V1, e$V2, type="l", main="98% clustering", xlab="% OTUs", ylab="Number of leaves", ylim=c(0,1600), col="red")
lines(f$V1, f$V2, type="l", col="blue")

#empty plot, placeholder for legend only
plot(c(0,1), c(0,1), type="n", axes=FALSE, main="", xlab="", ylab="")

legend("topleft", legend=c("100 bit score cutoff", "200 bit score cutoff"), col=c("red","blue"), bty="n", lty=1)

dev.off()

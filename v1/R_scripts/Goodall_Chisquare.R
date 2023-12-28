pdf(file="NMDS.pdf")
par(mfrow=c(2,2))

a<-read.table(file="Goodall.txt")
a
b<-as.matrix(a)
b
library(Matrix)
c<-tril(b)
c
d<-as.matrix(c)
d
e<-as.dist(d)
e

library(ecodist)
f<-nmds(e, nits=1000, epsilon=1e-12)
f
g<-nmds.min(f)
g
#pdf(file="Goodall.pdf")
plot(g, main="Goodall similarity metric", xlab="NMDS1", ylab="NMDS2", pch=c(15,15,15,15,17,17,17,17), col=c("steelblue","firebrick","olivedrab4", "mediumpurple4","steelblue","firebrick","olivedrab4", "mediumpurple4" ))

#print out lowest stress and highest r2
min(f$stress)
max(f$r2)

#dev.off()
#save.image(file="Goodall.RData")
#quit()

a<-read.table(file="ChiSquared.txt")
b<-as.matrix(a)
library(Matrix)
c<-tril(b)
d<-as.matrix(c)
e<-as.dist(d)
library(ecodist)
f<-nmds(e, nits=1000, epsilon=1e-12)
g<-nmds.min(f)
#pdf(file="Chi_squared.pdf")
plot(g, main="Chi squared", xlab="NMDS1", ylab="NMDS2", pch=c(15,15,15,15,17,17,17,17), col=c("steelblue","firebrick","olivedrab4", "mediumpurple4","steelblue","firebrick","olivedrab4", "mediumpurple4" ))

#print out lowest stress and highest r2
min(f$stress)
max(f$r2)

dev.off()
save.image(file="NMDS.RData")
quit()

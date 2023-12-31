pdf(file="NMDS.pdf")
par(mfrow=c(2,2))

a<-read.table(file="Goodall.txt")
b<-as.matrix(a)
library(Matrix)
c<-tril(b)
d<-as.matrix(c)
e<-as.dist(d)
library(ecodist)
f<-nmds(e, nits=1000)
g<-nmds.min(f)
#pdf(file="Goodall.pdf")
plot(g, main="Goodall similarity metric", xlab="NMDS1", ylab="NMDS2", pch=c(15,15,15,15,17,17,17,17), col=c("steelblue","firebrick","olivedrab4", "mediumpurple4","steelblue","firebrick","olivedrab4", "mediumpurple4" ))
#dev.off()
#save.image(file="Goodall.RData")
#quit()

A<-read.table(file="ChiSquared.txt")
B<-as.matrix(A)
library(Matrix)
C<-tril(B)
D<-as.matrix(C)
E<-as.dist(D)
library(ecodist)
F<-nmds(E, nits=1000)
G<-nmds.min(F)
#pdf(file="Chi_squared.pdf")
plot(G, main="Chi squared", xlab="NMDS1", ylab="NMDS2", pch=c(15,15,15,15,17,17,17,17), col=c("steelblue","firebrick","olivedrab4", "mediumpurple4","steelblue","firebrick","olivedrab4", "mediumpurple4" ))

dev.off()
save.image(file="NMDS.RData")
quit()

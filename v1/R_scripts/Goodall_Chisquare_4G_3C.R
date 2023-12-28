#pdf(file="NMDS.pdf")
#par(mfrow=c(2,2))

a<-read.table(file="Goodall.txt")
b<-as.matrix(a)
library(Matrix)
c<-tril(b)
d<-as.matrix(c)
e<-as.dist(d)
library(ecodist)
f<-nmds(e, maxdim=4)

f
min(f$stress)
max(f$r2)

g<-nmds.min(f, dims=0) #plot best overall
g

#plot scree
#stress<-f$stress
#plot(stress)
#axis.seq <- c(seq(1, 1, length = 10), seq(2, 2, length = 10),seq(3, 3, length = 10), seq(4, 4, length = 10), seq(5, 5,length = 10), seq(6, 6, length = 10), seq(7, 7, length = 10),seq(8, 8, length = 10), seq(9, 9, length = 10), seq(10, 10,length = 10))
#plot(stress ~ factor(axis.seq), main="Goodall scree")

pdf(file="Goodall.pdf")
#plot(g)
plot(g, main="Goodall similarity metric", pch=c(15,15,15,15,17,17,17,17), col=c("steelblue","firebrick","olivedrab4", "mediumpurple4","steelblue","firebrick","olivedrab4", "mediumpurple4" ))
dev.off()
#save.image(file="Goodall.RData")
#quit()

A<-read.table(file="ChiSquared.txt")
B<-as.matrix(A)
library(Matrix)
C<-tril(B)
D<-as.matrix(C)
E<-as.dist(D)
library(ecodist)
F<-nmds(E, maxdim=3)

F
min(F$stress)
max(F$r2)

G<-nmds.min(F, dims=0) #plot best overall
G

#plot scree
#stress2<-F$stress
#plot(stress2)
#axis.seq2 <- c(seq(1, 1, length = 10), seq(2, 2, length = 10),seq(3, 3, length = 10), seq(4, 4, length = 10), seq(5, 5,length = 10), seq(6, 6, length = 10), seq(7, 7, length = 10),seq(8, 8, length = 10), seq(9, 9, length = 10), seq(10, 10,length = 10))
#plot(stress2 ~ factor(axis.seq2), main="Chi squared scree")

pdf(file="Chi_squared.pdf")
#plot(G)
plot(G, main="Chi squared", pch=c(15,15,15,15,17,17,17,17), col=c("steelblue","firebrick","olivedrab4", "mediumpurple4","steelblue","firebrick","olivedrab4", "mediumpurple4" ))

dev.off()
save.image(file="NMDS.RData")
quit()

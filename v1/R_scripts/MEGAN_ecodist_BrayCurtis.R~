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

#adjust dims for scree
#f<-nmds(e, maxdim=10, nits=10)
#g<-nmds.min(f, dims=0) #best overall

#pdf(file="NMDS_scree.pdf")

#plot scree
#stress<-f$stress
#axis.seq <- c(seq(1, 1, length = 10), seq(2, 2, length = 10),seq(3, 3, length = 10), seq(4, 4, length = 10), seq(5, 5,length = 10), seq(6, 6, length = 10), seq(7, 7, length = 10),seq(8, 8, length = 10), seq(9, 9, length = 10), seq(10, 10,length = 10))
#plot(stress ~ factor(axis.seq), main="Scree")

#dev.off()

#do real NMDS now
h<-nmds(e, nits=100, maxdim=3)
i<-nmds.min(h)

pdf(file="NMDS.pdf")

library(RColorBrewer)
warm<-brewer.pal(9,"YlOrRd")
warm<-tail(warm,8)
warm
cool<-brewer.pal(9,"PuBuGn")
cool<-tail(cool,8)
cool

#pch 15 filled squares for rbcLF, pch 17 filled triangles for rbcLR
#col warm for PAD, col cool for WC
plot(i, main="Bray-Curtis disimilarity", pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17), col=c(warm,warm,cool,cool))

dev.off()

#create empty plot so legend prints separately
pdf(file="NMDS_legend.pdf")

plot(1,type="n",axes=FALSE,xlab="",ylab="")

legend("topleft", legend=labels, ncol=4, cex=0.75, pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17), col=c(warm,warm,cool,cool), bty="n")

save.image(file="NMDS.RData")
dev.off()
quit()

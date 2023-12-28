#read in sample by species matrix from MEGAN+Excel
FO<-read.csv("fungiorder.csv",head=TRUE,row.names=1)
FO
is.data.frame(FO)

#read in environment file
FO.env<-read.csv("sedimentfungi.env",head=TRUE,row.names=1)
FO.env
is.data.frame(FO.env)
attach(FO.env)
	
library(vegan)

#do NMDS for 2 dimensions
pdf("ordination2dim.pdf")
nmds2<-metaMDS(FO,k=2,trymax=100)
nmds2
plot(nmds2,type="n")
text(nmds2, display="sites",cex=0.8,col="red")
points(nmds2, display="spec", cex=0.7,pch=21,col="blue")

#fit environmental variables
nmds2.fit<-envfit(nmds2 ~ PAH + Alk + UCM, data=FO.env, perm=1000)
nmds2.fit
plot(nmds2.fit)
dev.off()

#do NMDS for 3 dimensions
pdf("ordination3dim.pdf")
par(mfrow=c(2,3))
nmds3<-metaMDS(FO,k=3,trymax=100)
nmds3

#plot NMDS1 vs NMDS2
plot(nmds3, type="n", choices=c(1,2))
text(nmds3, display="sites",cex=0.8,col="red", choices=c(1,2))
points(nmds3, display="spec", cex=0.7,pch=21,col="blue", choices=c(1,2))

#fit environmental variables
nmds3.fit12<-envfit(nmds3 ~ PAH + Alk + UCM, data=FO.env, perm=1000, choices=c(1,2))
nmds3.fit12
plot(nmds3.fit12)

#plot NMDS2 vs NMDS3
plot(nmds3, type="n", choices=c(2,3))
text(nmds3, display="sites",cex=0.8,col="red",choices=c(2,3))
points(nmds3, display="spec", cex=0.7,pch=21,col="blue",choices=c(2,3))

#fit environmental variables
nmds3.fit23<-envfit(nmds3 ~ PAH + Alk + UCM, data=FO.env, perm=1000,choices=c(2,3))
nmds3.fit23
plot(nmds3.fit23)

#plot NMDS3 vs NMDS1
plot(nmds3, type="n", choices=c(3,1))
text(nmds3, display="sites",cex=0.8,col="red",choices=c(3,1))
points(nmds3, display="spec", cex=0.7,pch=21,col="blue",choices=c(3,1))

#fit environmental variables
nmds3.fit31<-envfit(nmds3 ~ PAH + Alk + UCM, data=FO.env, perm=1000,choices=c(3,1))
nmds3.fit31
plot(nmds3.fit31)
dev.off()
	
#stressplot Shephards curve
pdf("stressplot2.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

pdf("stressplot3.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()

#directly compare a 2 dim versus 3 dim ordination
pdf("compare.pdf")
plot(procrustes(nmds2, nmds3))
dev.off()

#new function edited from http://www.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm
pdf("scree.pdf")
nmds.scree<-function(x) {
	plot(rep(1,10),replicate(10,metaMDS(x, autotransform=F,k=1)$stress/100),xlim=c(1,nrow(x)),ylim=c(0,0.5),xlab="Number dimensions",ylab="Stress",main="NMDS stress plot")

for(i in 1:(nrow(x)-2)) {
	points(rep(i+1,10),replicate(10,metaMDS(x, autotransform=F,k=i+1)$stress/100))
}
}

#make scree plot
FO
nmds.scree(FO)
dev.off()

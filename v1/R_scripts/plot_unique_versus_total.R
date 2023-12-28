a<-read.table("table.txt",header=TRUE,sep=",")
a
pdf("unique.pdf")
plot(a$uniques,a$total, main="Soil: 6 markers combined", xlab="OTUs unique to a site (%)", ylab="Total number of OTUS per site", pch=15, cex=0.6, col="black")
text(a$uniques,a$total,a$X, cex=0.6, pos=3, col=c("red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue"))
dev.off()

a<-read.table("table.txt",header=TRUE,sep=",")
a

pdf("uniques_all_norm.pdf")

par(mfrow=c(2,2))

library(RColorBrewer)
warm<-brewer.pal(9,"YlOrRd")
warm<-tail(warm,8)
warm<-rev(warm)
warm
cool<-brewer.pal(9,"PuBuGn")
cool<-tail(cool,8)
cool<-rev(cool)
cool

plot(a$total,a$uniques, main="Soil: 6 markers combined", xlab="% total OTUs per site", ylab="% unique per subsite", pch=22, cex=0.6, col="black", bg=c(warm,cool), xlim=c(0,100), ylim=c(0,100))

labels<-c("PAD1 (N=4,958)","PAD11 (N=4,656)","PAD14 (N=5,247)","PAD3 (N=4,682)","PAD33 (N=8,517)","PAD37 (N=5,216)","PAD38 (N=5,205)","PAD4 (N=6,301)","WC1 (N=5,689)","WC2 (N=5,714)","WC3 (N=6,014)","WC4 (N=5,931)","WC5 (N=5,260)","WC6 (N=5,994)","WC7 (N=5,590)","WC8 (N=5,269)")

plot(1,1,type="n",main="",xlab="", ylab="",axes=FALSE)

legend("topleft",legend=labels, cex=0.6, fill=c(warm,cool), bty="n")
dev.off()

#create pie charts
rbcL <- read.table(file="rbcL.txt",header=TRUE,sep=",")
rbcL
is.list(rbcL)

PAD_slices <- rbcL$PAD
PAD_slices

WC_slices <- rbcL$WC
WC_slices

lbls <- rbcL$X
lbls

PAD_pct <- round(PAD_slices/sum(PAD_slices)*100)
PAD_pct

PAD_lbls <- paste(lbls, PAD_pct, sep="\n")
PAD_lbls <- paste(PAD_lbls, "%", sep="")
PAD_lbls <- paste(PAD_lbls, "(", sep=" ")
PAD_lbls <- paste(PAD_lbls, PAD_slices, sep="")
PAD_lbls <- paste(PAD_lbls, ")", sep="")

WC_pct <- round(WC_slices/sum(WC_slices)*100)
WC_pct

WC_lbls <- paste(lbls, WC_pct, sep="\n")
WC_lbls <- paste(WC_lbls, "%", sep="")
WC_lbls <- paste(WC_lbls, "(", sep=" ")
WC_lbls <- paste(WC_lbls, WC_slices, sep="")
WC_lbls <- paste(WC_lbls, ")", sep="")

library(RColorBrewer)
#only need 4 colors this time
colors1 = brewer.pal(9,"Set1")
colors1
colors=colors1[1:4]
colors
#colors2 = brewer.pal(2,"Set2")
#colors = c(colors1, colors2)

pdf(file="rbcL_pie.pdf")
#par(mfrow=c(2,2), cex=0.5)
layout(matrix(c(1,2,3,3), 2,2, byrow=TRUE))

pie(PAD_slices, labels=PAD_lbls, main="PAD", cex.main=0.9, col=colors, clockwise=TRUE, radius=0.6, cex=0.5)
pie(WC_slices, labels=WC_lbls, main="WC", col=colors, cex.main=0.9, clockwise=TRUE, radius=0.6, cex=0.5)

#create bar chart
prop <-read.table(file="rbcL_taxa.txt", header=TRUE, sep=",", row.names=1)
prop
is.list(prop)

v<-cbind(prop$PAD, prop$WC)	
v

names<-row.names(prop)
names

library(RColorBrewer)
#need 11 colors
colors1 = brewer.pal(9,"Set1")
colors2 = brewer.pal(8,"Set2")
colors = c(colors1,colors2[1:2])
colors

par(xpd=T, mar=par()$mar+c(0,0,0,4))

barplot(v, beside=FALSE, names.arg=c("PAD (N=2348)", "WC (N=2829)"), cex.names=0.9, ylab="OTUs (%)", col=colors, width=0.5, ylim=c(0,100), xlim=c(0,1.25), space=0.1)

legend(1.2,100, legend=rev(names), cex=0.5, fill=rev(colors), bty="n")

par(mar=c(5,4,4,2)+0.1)

dev.off()

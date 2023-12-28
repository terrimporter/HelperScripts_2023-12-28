abund <- read.table(file="abund.csv", header=TRUE, sep=",")
abund
is.list(abund)

v <- cbind(abund$GenBank, abund$GenBank_barcode, abund$GenBank_family)
v

library(RColorBrewer)
#need 29 colors!
colors1 = brewer.pal(9,"Set1")
colors2 = brewer.pal(8,"Set2")
colors3 = brewer.pal(12,"Set3")
colors = c(colors1, colors2, colors3)

pdf(file="stacked.pdf")

par(xpd=T, mar=par()$mar+c(0,0,0,4))

barplot(v, beside=FALSE, names.arg=c("GenBank", "GenBank-barcode", "GenBank-family"), cex.names=0.5, ylim=c(0,300000), ylab="Number of sequences", col=colors, width = 1, xlim = c(0,6), space=0.1)

legend(4,300000,legend=rev(abund$X), cex=0.5, fill=rev(colors))

par(mar=c(5,4,4,2)+0.1)
dev.off()

parta <- read.table(file="4a.csv", header=TRUE, sep=",")
parta
is.list(parta)

partb <- read.table(file="4b.csv", header=TRUE, sep=",")
partb

v <- cbind(parta$Morphology, parta$GenBank_NBC_MEGAN, parta$GenBank_barcode_NBC_MEGAN, parta$GenBank_family_NBC_MEGAN, parta$BLAST_MEGAN)
v

w <- cbind(partb$Morphology, partb$GenBank_NBC_MEGAN_20, partb$GenBank_barcode_NBC_MEGAN_0, partb$GenBank_family_NBC_MEGAN_60)
w

library(RColorBrewer)
#need 21 colors!
colors1 = brewer.pal(9,"Set1")
colors2 = brewer.pal(8,"Set2")
colors3 = brewer.pal(4,"Set3")
colors = c(colors1, colors2, colors3)

pdf(file="stacked.pdf")

par(mfcol=c(2,1))

par(xpd=T, mar=par()$mar+c(0,0,0,4))

barplot(v, beside=FALSE, names.arg=c("Morphology","GenBank", "GenBank-barcode", "GenBank-family", "BLAST"), cex.names=0.5, ylim=c(0,1000), ylab="Number of sequences", col=colors, width = 1, xlim = c(0,6), space=0.1)

legend(6,1000,legend=rev(parta$X), cex=0.5, fill=rev(colors))

barplot(w, beside=FALSE, names.arg=c("Morphology","GenBank\n(20%)", "GenBank-barcode\n(0%)", "GenBank-family\n(60%)"), cex.names=0.5, ylim=c(0,1000), ylab="Number of sequences", col=colors, width = 1, xlim=c(0,6), space=0.1)

par(mar=c(5,4,4,2)+0.1)
dev.off()

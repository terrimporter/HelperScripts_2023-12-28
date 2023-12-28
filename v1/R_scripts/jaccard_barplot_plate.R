a<-read.table(file="a.txt",header=TRUE, sep="\t", row.names=1)
a
b<-read.table(file="b.txt",header=TRUE, sep="\t", row.names=1)
b
c<-read.table(file="c.txt",header=TRUE, sep="\t", row.names=1)
c
d<-read.table(file="d.txt",header=TRUE, sep="\t", row.names=1)
d
e<-read.table(file="e.txt",header=TRUE, sep="\t", row.names=1)
e
f<-read.table(file="f.txt",header=TRUE, sep="\t", row.names=1)
f

A<-cbind(a$Within, a$Among)
A
At<-t(A) #transpose

Astdev<-cbind(a$WithinSTDEV, a$AmongSTDEV)
Astdevt<-t(Astdev) #transpose

Atogether<-At+Astdevt
Atogether #height of bars (arrows)

B<-cbind(b$Within, b$Among)
B
Bt<-t(B)

Bstdev<-cbind(b$WithinSTDEV, b$AmongSTDEV)
Bstdevt<-t(Bstdev)

Btogether<-Bt+Bstdevt

C<-cbind(c$Within, c$Among)
C
Ct<-t(C)

Cstdev<-cbind(c$WithinSTDEV, c$AmongSTDEV)
Cstdevt<-t(Cstdev)

Ctogether<-Ct+Cstdevt

D<-cbind(d$X0_pooled, d$X3_pooled, d$X9_pooled)
D
Dt<-t(D)

E<-cbind(e$X0_pooled, e$X3_pooled, e$X9_pooled)
E
Et<-t(E)

F<-cbind(f$X0_pooled, f$X3_pooled, f$X9_pooled)
F
Ft<-t(F)

library(RColorBrewer)
colors1 = brewer.pal(9, "Set1")
col_within_among = colors1[1:2]
col_within_among
col_0_3_9_pooled = colors1[3:5]
col_0_3_9_pooled

pdf(file="jaccard_barplot_plate.pdf")
par(mfrow=c(3,3))
#par(xpd=T, mar=par()$mar+c(0,0,0,4))

plotA<-barplot(as.matrix(At), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="No pooling\n(N=144)", col=col_within_among, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotA, as.matrix(At), plotA, as.matrix(Atogether), angle=90, code=2, length=0.02)

plotB<-barplot(as.matrix(Bt), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="3 pooled cores per sample\n(N=48)", col=col_within_among, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotB, as.matrix(Bt), plotB, as.matrix(Btogether), angle=90, code=2, length=0.02)

plotC<-barplot(as.matrix(Ct), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="9 pooled cores per sample\n(N=16)", col=col_within_among, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotC, as.matrix(Ct), plotC, as.matrix(Ctogether), angle=90, code=2, length=0.02)

plotD<-barplot(as.matrix(Dt), beside=TRUE, ylab="Average Jaccard dissimilarity difference\n(Among - within sites)", ylim=c(0,0.16), cex.lab=0.5, main="Jaccard dissimilarity difference", col=col_0_3_9_pooled, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

plotE<-barplot(as.matrix(Et), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.0), cex.lab=0.9, main="Jaccard dissimilarity\nWithin sites", col=col_0_3_9_pooled, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

plotF<-barplot(as.matrix(Ft), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.0), cex.lab=0.9, main="Jaccard dissimilarity\nAmong sites", col=col_0_3_9_pooled, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")

legend(0, 1, legend=c("Within sites", "Among sites", "No pooling", "3 cores pooled per sample", "9 cores pooled per sample"), cex=0.9, fill=c(col_within_among,col_0_3_9_pooled), bty="n")

#par(mar=c(5,4,4,2)+0.1)
dev.off()

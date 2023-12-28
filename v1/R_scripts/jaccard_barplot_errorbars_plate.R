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

Estdev<-cbind(e$X0_pooledSTDEV, e$X3_pooledSTDEV, e$X9_pooledSTDEV)
Estdevt<-t(Estdev)

Etogether<-Et+Estdevt

F<-cbind(f$X0_pooled, f$X3_pooled, f$X9_pooled)
F
Ft<-t(F)

Fstdev<-cbind(f$X0_pooledSTDEV, f$X3_pooledSTDEV, f$X9_pooledSTDEV)
Fstdevt<-t(Fstdev)

Ftogether<-Ft+Fstdevt

library(RColorBrewer)
colors1 = brewer.pal(9, "Set1")
col_within_among = colors1[1:2]
col_within_among
col_0_3_9_pooled = colors1[3:5]
col_0_3_9_pooled

pdf(file="jaccard_barplot_plate.pdf")
par(mfrow=c(3,3))
#par(xpd=T, mar=par()$mar+c(0,0,0,4))

plotA<-barplot(as.matrix(At), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="No pooling", col=col_within_among, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotA, as.matrix(At), plotA, as.matrix(Atogether), angle=90, code=2, length=0.02)

plotB<-barplot(as.matrix(Bt), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="3 pooled soil samples per block", col=col_within_among, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotB, as.matrix(Bt), plotB, as.matrix(Btogether), angle=90, code=2, length=0.02)

plotC<-barplot(as.matrix(Ct), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="3 pooled blocks per subsite", col=col_within_among, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotC, as.matrix(Ct), plotC, as.matrix(Ctogether), angle=90, code=2, length=0.02)

plotD<-barplot(as.matrix(Dt), beside=TRUE, ylab="Average Jaccard dissimilarity difference\n(Among - within sites)", ylim=c(0,0.16), cex.lab=0.5, main="Jaccard dissimilarity difference", col=col_0_3_9_pooled, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

plotE<-barplot(as.matrix(Et), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="Jaccard dissimilarity\nWithin sites", col=col_0_3_9_pooled, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotE, as.matrix(Et), plotE, as.matrix(Etogether), angle=90, code=2, length=0.02)

plotF<-barplot(as.matrix(Ft), beside=TRUE, ylab="Jaccard dissimilarity", ylim=c(0,1.2), cex.lab=0.9, main="Jaccard dissimilarity\nAmong sites", col=col_0_3_9_pooled, names.arg=c("rbcLF", "rbcLR", "ITS1", "ITS4", "16V3", "16V6"), cex.names=0.9, las=2, cex.axis=0.9)

arrows(plotF, as.matrix(Ft), plotF, as.matrix(Ftogether), angle=90, code=2, length=0.02)

plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")

legend(0, 1, legend=c("Within PAD or WC sites", "Among PAD and WC sites", "No pooling", "3 pooled soil samples per block", "3 pooled blocks per subsite"), cex=0.8, fill=c(col_within_among,col_0_3_9_pooled), bty="n")

#par(mar=c(5,4,4,2)+0.1)
dev.off()

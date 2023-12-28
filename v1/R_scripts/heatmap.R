a<-read.table("heatmap.txt", header=T)
a
b<-as.matrix(a,nrow= 305,ncol= 24)
b
c<-t(b)
c
f<-c[,order(colSums(c),decreasing=TRUE)]
f
g<-f[(order(rowSums(f),decreasing=TRUE)), ]
g

is.matrix(g)
is.numeric(g)

library(RColorBrewer)
golden<-brewer.pal(5,"YlOrRd")

library(gplots)
pdf()
heatmap.2(g,Rowv=NA, Colv=NA, col=golden, xlab="Orthologs", ylab="Taxa", Key=T, Keysize=1.5, density.info="none",trace="none")
dev.off()


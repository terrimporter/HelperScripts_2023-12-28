a<-read.table("phylum_cmp.megan-chart", sep="\t", header=TRUE)
a
b<-as.matrix(a[,-1])
b
c<-data.matrix(b)
c
#d<-c[,-1] #exclude first column w/ labels
rowNames<-a[,1] #keep row names

#c<-t(b)
#c
#f<-c[,order(colSums(c),decreasing=TRUE)]
#f
#g<-f[(order(rowSums(f),decreasing=TRUE)), ]
#g

#is.matrix(c)
#is.numeric(c)
#d<-data.matrix(c,rownames.force=NA)
#d
is.matrix(c)
is.numeric(c)

library(RColorBrewer)
reds<-brewer.pal(9,"Reds")

library(gplots)
pdf()
heatmap.2(c,Rowv=NA, Colv=TRUE, distfun=dist, hclustfun=hclust, dendrogram="column", symm=FALSE, col=reds, xlab="Subsites", ylab="Taxa (Phylum)", Key=T, Keysize=0.25, density.info="none",trace="none", labRow=rowNames)
dev.off()


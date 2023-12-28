##### Read infiles #####

#Jaccard pairwise comparisons among blocks
blocks <- read.table(file="jaccard_blocks.txt", header=FALSE)
blocks
is.list(blocks)

##### Frequency histograms #####

pdf(file="jaccard_hist.pdf")

par(mfcol=c(2,2))

#GenBank classifier
a<-hist(blocks$V1, breaks=seq(0,1.0,by=0.1), main="Soil, Subsampled", xlab="Jaccard dissimilarity", ylab="Soil block pairs")
a$density <- (a$counts/sum(a$counts))*100
a

#dev.off()

##### Percent histograms #####

#pdf(file="jaccard_percent_histogram.pdf")

#GenBank classifier
plot(a, freq=FALSE, main="Soil, Subsampled", xlab="Jaccard dissimilarity", ylab="Soil block pairs (%)", ylim=c(0,100))

dev.off()

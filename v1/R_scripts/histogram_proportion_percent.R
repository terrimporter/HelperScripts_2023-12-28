##### Read infiles #####

#GenBank classifier
c100_v2 <- read.table(file="C100_v2_family.txt", header=FALSE)
c100_v2
is.list(c100_v2)

c200_v2 <- read.table(file="C200_v2_family.txt", header=FALSE)
c200_v2

Reference_v2 <- read.table(file="Reference_v2_family.txt",header=FALSE)
Reference_v2

#GenBank-barcode classifier
c100_v3 <- read.table(file="C100_v3_family.txt", header=FALSE)
c100_v3

c200_v3 <- read.table(file="C200_v3_family.txt", header=FALSE)
c200_v3

Reference_v3 <- read.table(file="Reference_v3_family.txt",header=FALSE)
Reference_v3

#GenBank-family classifier
c100_v4 <- read.table(file="C100_v4_family.txt", header=FALSE)
c100_v4

c200_v4 <- read.table(file="C200_v4_family.txt", header=FALSE)
c200_v4

Reference_v4 <- read.table(file="Reference_v4_family.txt",header=FALSE)
Reference_v4

##### Frequency histograms #####

pdf(file="C_freq_histogram.pdf")

par(mfcol=c(3,3))

#GenBank classifier
a<-hist(c100_v2$V1, breaks=seq(0,1.0,by=0.1), main="C100_v2 (N=9880)", xlab="Bootstrap support", ylab="Frequency")
a$density <- (a$counts/sum(a$counts))*100
a

b<-hist(c200_v2$V1, breaks=seq(0,1.0,by=0.1), main="C200_v2 (N=9323)", xlab="Bootstrap support", ylab="Frequency")
b$density<-(b$counts/sum(b$counts))*100
b

c<-hist(Reference_v2$V1, breaks=seq(0,1.0,by=0.1), main="Reference_v2 (N=286)", xlab="Bootstrap support", ylab="Frequency")
c$density<-(c$counts/sum(c$counts))*100
c

#GenBank-barcode classifier
d<-hist(c100_v3$V1, breaks=seq(0,1.0,by=0.1), main="C100_v3", xlab="Bootstrap support", ylab="Frequency")
d$density<-(d$counts/sum(d$counts))*100

e<-hist(c200_v3$V1, breaks=seq(0,1.0,by=0.1), main="C200_v3", xlab="Bootstrap support", ylab="Frequency")
e$density<-(e$counts/sum(e$counts))*100

f<-hist(Reference_v3$V1, breaks=seq(0,1.0,by=0.1), main="Reference_v3", xlab="Bootstrap support", ylab="Frequency")
f$density<-(f$counts/sum(f$counts))*100

#GenBank-family classifier
g<-hist(c100_v4$V1, breaks=seq(0,1.0,by=0.1), main="C100_v4", xlab="Bootstrap support", ylab="Frequency")
g$density<-(g$counts/sum(g$counts))*100

h<-hist(c200_v4$V1, breaks=seq(0,1.0,by=0.1), main="C200_v4", xlab="Bootstrap support", ylab="Frequency")
h$density<-(h$counts/sum(h$counts))*100

i<-hist(Reference_v4$V1, breaks=seq(0,1.0,by=0.1), main="Reference_v4", xlab="Bootstrap support", ylab="Frequency")
i$density<-(i$counts/sum(i$counts))*100

dev.off()

##### Percent histograms #####

pdf(file="C_percent_histogram.pdf")

par(mfcol=c(3,3))

#GenBank classifier
plot(a, freq=FALSE, main="C100_v2 (N=9692, average length = 155 bp)", cex.main = 0.5, xlab="", ylab="Sequences (%)", ylim=c(0,100))

plot(b, freq=FALSE, main="C200_v2 (N=9275, average length = 254 bp)", cex.main = 0.5, xlab="", ylab="Sequences (%)", ylim=c(0,100))

plot(c, freq=FALSE, main="Reference_v2 (N=282, average length = 626 bp)", cex.main = 0.5, xlab="Bootstrap support", ylab="Sequences (%)", ylim=c(0,100))

#GenBank-barcode classifier
plot(d, freq=FALSE, main="C100_v3", cex.main=0.5, xlab="", ylab="", ylim=c(0,100))

plot(e, freq=FALSE, main="C200_v3", cex.main=0.5, xlab="", ylab="", ylim=c(0,100))

plot(f, freq=FALSE, main="Reference_v3", cex.main=0.5, xlab="Bootstrap support", ylab="", ylim=c(0,100))

#GenBank-family classifier
plot(g, freq=FALSE, main="C100_v4", cex.main=0.5, xlab="", ylab="", ylim=c(0,100))

plot(h, freq=FALSE, main="C200_v4", cex.main=0.5, xlab="", ylab="", ylim=c(0,100))

plot(i, freq=FALSE, main="Reference_v4", cex.main=0.5, xlab="Bootstrap support",ylab="", ylim=c(0,100))

dev.off()

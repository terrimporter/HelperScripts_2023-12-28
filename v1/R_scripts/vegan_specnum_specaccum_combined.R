ITS4F<-read.csv(file="incidence.csv", head=TRUE, row.names=1)
ITS4F

ITS4F.env<-read.table(file="ITS4F.env",head=TRUE)
ITS4F.env

library(vegan)

#specnum<-specnumber(ITS4F)
#specnum

#pdf("specnum_boxplot.pdf")

#boxplot(specnum ~ ITS4F.env$Site, ylab="Number of OTUs")

#dev.off()

#t.test(specnum ~ ITS4F.env$Site)

#pdf("specaccum_samplingcurve.pdf")

#par(mfrow=c(2,2))

#edit this using marker.env file in excel
ITS4F_PAD<-ITS4F[1:69,]
ITS4F_WC<-ITS4F[70:135,]

#remove columns with only zeros
ITS4F_PAD<-ITS4F_PAD[,colSums(ITS4F_PAD) !=0]
ITS4F_WC<-ITS4F_WC[,colSums(ITS4F_WC) !=0]

#plot(specaccum(ITS4F_PAD), xlab="Number of soil samples", ylab="Number of OTUs", main="PAD")
#plot(specaccum(ITS4F_WC), xlab="Number of soil samples", ylab="Number of OTUs", main="WC")

#dev.off()

pdf("ITS4F_combined_sampling_curve.pdf")

plot(specaccum(ITS4F_PAD), xlab="Number of soil samples", ylab="Number of OTUs", main="ITS4F", col="red", xlim=c(0,70), ylim=c(0,4000))
par(new=T)
plot(specaccum(ITS4F_WC), xlab="", ylab="", main="", col="blue", xlim=c(0,70), ylim=c(0,4000))

legend("bottomright", legend=c("PAD","WC"), col=c("red","blue"), pch="-", bty="n")

dev.off()

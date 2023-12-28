misclass0 <- read.table(file="misclass0.csv",header=TRUE,sep=",")
misclass0

misclass0_Raphidioptera<-misclass0[24,]
misclass0_Mantodea<-misclass0[13,]
misclass0_Mecoptera<-misclass0[15,]
misclass0_Psocoptera<-misclass0[23,]
misclass0_undef_Insecta<-misclass0[29,]
misclass0_Neuroptera<-misclass0[17,]
misclass0_Siphonaptera<-misclass0[25,]
misclass0_Blattodea<-misclass0[2,]

misclass0<-misclass0[-c(2,13,15,17,23,24,25,29),]

is.list(misclass0)

misclass90 <- read.table(file="misclass90.csv",header=TRUE,sep=",")
misclass90

misclass90_Raphidioptera<-misclass90[24,]
misclass90_Mantodea<-misclass90[13,]
misclass90_Mecoptera<-misclass90[15,]
misclass90_Psocoptera<-misclass90[23,]
misclass90_undef_Insecta<-misclass90[29,]
misclass90_Neuroptera<-misclass90[17,]
misclass90_Siphonaptera<-misclass90[25,]
misclass90_Blattodea<-misclass90[2,]

misclass90<-misclass90[-c(2,13,15,17,23,24,25,29),]


pdf(file="Misclassified.pdf")

par(mfrow=c(2,1))
par(xpd=T, mar=par()$mar+c(0,0,0,4))

library(RColorBrewer)
colors1=brewer.pal(8,"Set1")

#Plot misclassified, no cutoff
plot(x=misclass0$Total, y=misclass0$Misclassified, log="x", axes=FALSE, ylim=c(0,100), xlab="Log number of sequences", ylab="% Misclassified (No bootstrap support cutoff)", pch=20)
lines(x=misclass0_Raphidioptera$Total,y=misclass0_Raphidioptera$Misclassified,type="p", pch=20, col=colors1[1])
lines(x=misclass0_Mantodea$Total,y=misclass0_Mantodea$Misclassified,type="p", pch=20, col=colors1[2])
lines(x=misclass0_Mecoptera$Total,y=misclass0_Mecoptera$Misclassified,type="p", pch=20, col=colors1[3])
lines(x=misclass0_Psocoptera$Total,y=misclass0_Psocoptera$Misclassified,type="p", pch=20, col=colors1[4])
lines(x=misclass0_undef_Insecta$Total,y=misclass0_undef_Insecta$Misclassified,type="p", pch=20, col=colors1[5])
lines(x=misclass0_Neuroptera$Total,y=misclass0_Neuroptera$Misclassified,type="p", pch=20, col=colors1[6])
lines(x=misclass0_Siphonaptera$Total,y=misclass0_Siphonaptera$Misclassified,type="p", pch=20, col=colors1[7])
lines(x=misclass0_Blattodea$Total,y=misclass0_Blattodea$Misclassified,type="p", pch=20, col=colors1[8])

#x-axis
xticks<-c(1,10,100,1000,10000,100000)
axis(1, at=xticks,labels=xticks,pos=0)

#y-axis
#need two axis commands, one for the line, one for the ticks and labels
axis(2,at=c(0,100), labels=c("",""),lwd.ticks=0, pos=1)
axis(2, at=seq(0,100,by=20), lwd=0, lwd.ticks=1, pos=1, cex.axis=0.75)

lines(x=c(1,100000), y=c(10,10), lty=2, col="black")

legend(100000, 100, legend=c("Raphidioptera, n=5", "Mantodea, n=74", "Mecoptera, n=34", "Psocopera, n=7", "undef_Insecta, n=9", "Neuroptera, n=454", "Siphonaptera, n=9", "Blattodea, n=173"), cex=0.5, fill=colors1)

#Plot misclassified, 90% cutoff
plot(x=misclass90$Total, y=misclass90$Misclassified, log="x", axes=FALSE, ylim=c(0,100), xlab="Log number of sequences", ylab="% Misclassified (90% bootstrap support cutoff", pch=20)
lines(x=misclass90_Raphidioptera$Total,y=misclass90_Raphidioptera$Misclassified,type="p", pch=20, col=colors1[1])
lines(x=misclass90_Mantodea$Total,y=misclass90_Mantodea$Misclassified,type="p", pch=20, col=colors1[2])
lines(x=misclass90_Mecoptera$Total,y=misclass90_Mecoptera$Misclassified,type="p", pch=20, col=colors1[3])
lines(x=misclass90_Psocoptera$Total,y=misclass90_Psocoptera$Misclassified,type="p", pch=20, col=colors1[4])
lines(x=misclass90_undef_Insecta$Total,y=misclass90_undef_Insecta$Misclassified,type="p", pch=20, col=colors1[5])
lines(x=misclass90_Neuroptera$Total,y=misclass90_Neuroptera$Misclassified,type="p", pch=20, col=colors1[6])
lines(x=misclass90_Siphonaptera$Total,y=misclass90_Siphonaptera$Misclassified,type="p", pch=20, col=colors1[7])
lines(x=misclass90_Blattodea$Total,y=misclass90_Blattodea$Misclassified,type="p", pch=20, col=colors1[8])

#x-axis
xticks<-c(1,10,100,1000,10000,100000)
axis(1, at=xticks,labels=xticks,pos=0)

#y-axis
#need two axis commands, one for the line, one for the ticks and labels
axis(2,at=c(0,100), labels=c("",""),lwd.ticks=0, pos=1)
axis(2, at=seq(0,100,by=20), lwd=0, lwd.ticks=1, pos=1, cex.axis=0.75)

lines(x=c(1,100000), y=c(10,10), lty=2, col="black")

par(mar=c(5,4,4,2)+0.1)

dev.off()

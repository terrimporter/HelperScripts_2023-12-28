time <- read.table(file="time.csv", header=TRUE, sep=",")
time
is.list(time)

#v<-cbind(correct$FULL, correct$X400, correct$X200, correct$X100, correct$X50)
#v

#library(RColorBrewer)
#only need 4 colors this time
#colors1 = brewer.pal(4,"Set1")
#colors2 = brewer.pal(8,"Set2")
#colors3 = brewer.pal(12,"Set3")
#colors = c(colors1, colors2, colors3)

pdf(file="time.pdf")

par(xpd=T, mar=par()$mar+c(0,0,0,4))

barplot(time$Classifications_minute, ylab="Queries processed per minute", ylim=c(0,5000), col="black", names.arg=c("BLAST", "NBC GenBank", "NBC GenBank-\nbarcode", "NBC GenBank-\nfamily"), cex.names=0.75, space=0.1)

#barplot(time, beside=TRUE, ylab="Correct (%)", yaxp=c(0,100,10), ylim=c(0,100), col=colors1, names.arg=c("500 bp+", "400 bp", "200 bp", "100 bp", "50 bp"))

#lines(x=c(0,25), y=c(95,95), lty=2)

#legend(25,100,legend=correct$X, cex=1, fill=colors1)

par(mar=c(5,4,4,2)+0.1)
dev.off()

misclass0 <- read.table(file="misclass0.csv",header=TRUE,sep=",")
misclass0

misclass0_Raphidioptera<-misclass0[1,]
misclass0_Craterostigmomorpha<-misclass0[2,]
misclass0_Siphonaptera<-misclass0[3,]
misclass0_Pauropoda<-misclass0[4,]
misclass0_Embioptera<-misclass0[5,]
misclass0_Hemiptera<-misclass0[6,]
misclass0_Geophilomorpha<-misclass0[7,]
misclass0_Neuroptera<-misclass0[9,]
misclass0_Thysanoptera<-misclass0[10,]
misclass0_Mantodea<-misclass0[11,]
misclass0_Scolopendromorpha<-misclass0[13,]
misclass0_Pantopoda<-misclass0[29,]
misclass0_Amblypygi<-misclass0[40,]
misclass0_Astigmata<-misclass0[43,]
misclass0_Blattodea<-misclass0[45,]
misclass0_Cumacea<-misclass0[49,]
misclass0_Dermaptera<-misclass0[50,]
misclass0_Endeostigmata<-misclass0[51,]
misclass0_Glomerida<-misclass0[53,]
misclass0_Grylloblattodea<-misclass0[54,]
misclass0_Halocryprida<-misclass0[55,]
misclass0_Harpacticoida<-misclass0[56,]
misclass0_Julida<-misclass0[58,]
misclass0_Lithobiomorpha<-misclass0[61,]
misclass0_Myodocopida<-misclass0[66,]
misclass0_Poecilostomatoida<-misclass0[73,]
misclass0_Polyzoniida<-misclass0[75,]
misclass0_Pseudoscorpiones<-misclass0[78,]
misclass0_Psocoptera<-misclass0[79,]
misclass0_Siphonophorida<-misclass0[82,]
misclass0_Spirostreptida<-misclass0[84,]
misclass0_Strepsiptera<-misclass0[85,]
misclass0_Tanaidacea<-misclass0[87,]
misclass0_undef_Diplopoda<-misclass0[89,]
misclass0_undef_Malacostraca<-misclass0[92,]

misclass0<-misclass0[-c(1,2,3,4,5,6,7,9,10,11,13,29,40,43,45,49,50,51,53,54,55,56,58,61,66,73,75,78,79,82,84,85,87,89,92),]

is.list(misclass0)

misclass90 <- read.table(file="misclass90.csv",header=TRUE,sep=",")
misclass90

misclass90_Raphidioptera<-misclass90[1,]
misclass90_Craterostigmomorpha<-misclass90[2,]
misclass90_Siphonaptera<-misclass90[3,]
misclass90_Pauropoda<-misclass90[4,]
misclass90_Embioptera<-misclass90[5,]
misclass90_Hemiptera<-misclass90[6,]
misclass90_Geophilomorpha<-misclass90[7,]
misclass90_Neuroptera<-misclass90[9,]
misclass90_Thysanoptera<-misclass90[10,]
misclass90_Mantodea<-misclass90[11,]
misclass90_Scolopendromorpha<-misclass90[13,]
misclass90_Pantopoda<-misclass90[29,]
misclass90_Amblypygi<-misclass90[40,]
misclass90_Astigmata<-misclass90[43,]
misclass90_Blattodea<-misclass90[45,]
misclass90_Cumacea<-misclass90[49,]
misclass90_Dermaptera<-misclass90[50,]
misclass90_Endeostigmata<-misclass90[51,]
misclass90_Glomerida<-misclass90[53,]
misclass90_Grylloblattodea<-misclass90[54,]
misclass90_Halocryprida<-misclass90[55,]
misclass90_Harpacticoida<-misclass90[56,]
misclass90_Julida<-misclass90[58,]
misclass90_Lithobiomorpha<-misclass90[61,]
misclass90_Myodocopida<-misclass90[66,]
misclass90_Poecilostomatoida<-misclass90[73,]
misclass90_Polyzoniida<-misclass90[75,]
misclass90_Pseudoscorpiones<-misclass90[78,]
misclass90_Psocoptera<-misclass90[79,]
misclass90_Siphonophorida<-misclass90[82,]
misclass90_Spirostreptida<-misclass90[84,]
misclass90_Strepsiptera<-misclass90[85,]
misclass90_Tanaidacea<-misclass90[87,]
misclass90_undef_Diplopoda<-misclass90[89,]
misclass90_undef_Malacostraca<-misclass90[92,]

misclass90<-misclass90[-c(1,2,3,4,5,6,7,9,10,11,13,29,40,43,45,49,50,51,53,54,55,56,58,61,66,73,75,78,79,82,84,85,87,89,92),]

pdf(file="Misclassified.pdf")

par(mfrow=c(2,1))

library(RColorBrewer)
colors1=brewer.pal(8,"Set1")
colors2=brewer.pal(8,"Set2")
colors3=brewer.pal(12,"Set3")
colors4=brewer.pal(7,"Dark2")
colors=c(colors1,colors2,colors3,colors4)

#Plot misclassified, no cutoff
plot(x=misclass0$Total, y=misclass0$Misclassified, log="x", axes=FALSE, ylim=c(0,100), xlim=c(1,100000), xlab="Number of sequences", ylab="% Misclassified", pch=21, cex=1.25, col="black", bg="black")
lines(x=misclass0_Raphidioptera$Total,y=misclass0_Raphidioptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[1])
lines(x=misclass0_Craterostigmomorpha$Total,y=misclass0_Craterostigmomorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[2])
lines(x=misclass0_Siphonaptera$Total,y=misclass0_Siphonaptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[3])
lines(x=misclass0_Pauropoda$Total,y=misclass0_Pauropoda$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[4])
lines(x=misclass0_Embioptera$Total,y=misclass0_Embioptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[5])
lines(x=misclass0_Hemiptera$Total,y=misclass0_Hemiptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[6])
lines(x=misclass0_Geophilomorpha$Total,y=misclass0_Geophilomorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[7])
lines(x=misclass0_Neuroptera$Total,y=misclass0_Neuroptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[8])
lines(x=misclass0_Thysanoptera$Total,y=misclass0_Thysanoptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[9])
lines(x=misclass0_Mantodea$Total,y=misclass0_Mantodea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[10])
lines(x=misclass0_Scolopendromorpha$Total,y=misclass0_Scolopendromorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[11])
lines(x=misclass0_Pantopoda$Total,y=misclass0_Pantopoda$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[12])
lines(x=misclass0_Amblypygi$Total,y=misclass0_Amblypygi$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[13])
lines(x=misclass0_Astigmata$Total,y=misclass0_Astigmata$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[14])
lines(x=misclass0_Blattodea$Total,y=misclass0_Blattodea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[15])
lines(x=misclass0_Cumacea$Total,y=misclass0_Cumacea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[16])
lines(x=misclass0_Dermaptera$Total,y=misclass0_Dermaptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[17])
lines(x=misclass0_Endeostigmata$Total,y=misclass0_Endeostigmata$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[18])
lines(x=misclass0_Glomerida$Total,y=misclass0_Glomerida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[19])
lines(x=misclass0_Grylloblattodea$Total,y=misclass0_Grylloblattodea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[20])
lines(x=misclass0_Halocryprida$Total,y=misclass0_Halocryprida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[21])
lines(x=misclass0_Harpacticoida$Total,y=misclass0_Harpacticoida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[22])
lines(x=misclass0_Julida$Total,y=misclass0_Julida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[23])
lines(x=misclass0_Lithobiomorpha$Total,y=misclass0_Lithobiomorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[24])
lines(x=misclass0_Myodocopida$Total,y=misclass0_Myodocopida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[25])
lines(x=misclass0_Poecilostomatoida$Total,y=misclass0_Poecilostomatoida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[26])
lines(x=misclass0_Polyzoniida$Total,y=misclass0_Polyzoniida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[27])
lines(x=misclass0_Pseudoscorpiones$Total,y=misclass0_Pseudoscorpiones$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[28])
lines(x=misclass0_Psocoptera$Total,y=misclass0_Psocoptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[29])
lines(x=misclass0_Siphonophorida$Total,y=misclass0_Siphonophorida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[30])
lines(x=misclass0_Spirostreptida$Total,y=misclass0_Spirostreptida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[31])
lines(x=misclass0_Strepsiptera$Total,y=misclass0_Strepsiptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[32])
lines(x=misclass0_Tanaidacea$Total,y=misclass0_Tanaidacea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[33])
lines(x=misclass0_undef_Diplopoda$Total,y=misclass0_undef_Diplopoda$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[34])
lines(x=misclass0_undef_Malacostraca$Total,y=misclass0_undef_Malacostraca$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[35])

#x-axis
xticks<-c(1,10,100,1000,10000,100000)
axis(1, at=xticks,labels=xticks,pos=0)

#y-axis
#need two axis commands, one for the line, one for the ticks and labels
axis(2, at=c(0,100), labels=c("",""),lwd.ticks=0, pos=1)
axis(2, at=seq(0,100,by=20), lwd=0, lwd.ticks=1, pos=1, cex.axis=0.75)

lines(x=c(1,100000), y=c(10,10), lty=2, col="black")

#Plot misclassified, 90% cutoff
plot(x=misclass90$Total, y=misclass90$Misclassified, log="x", axes=FALSE, ylim=c(0,100), xlim=c(1,100000), xlab="Number of sequences", ylab="% Misclassified", pch=21, cex=1.25, col="black", bg="black")
lines(x=misclass90_Raphidioptera$Total,y=misclass90_Raphidioptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[1])
lines(x=misclass90_Craterostigmomorpha$Total,y=misclass90_Craterostigmomorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[2])
lines(x=misclass90_Siphonaptera$Total,y=misclass90_Siphonaptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[3])
lines(x=misclass90_Pauropoda$Total,y=misclass90_Pauropoda$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[4])
lines(x=misclass90_Embioptera$Total,y=misclass90_Embioptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[5])
lines(x=misclass90_Hemiptera$Total,y=misclass90_Hemiptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[6])
lines(x=misclass90_Geophilomorpha$Total,y=misclass90_Geophilomorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[7])
lines(x=misclass90_Neuroptera$Total,y=misclass90_Neuroptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[8])
lines(x=misclass90_Thysanoptera$Total,y=misclass90_Thysanoptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[9])
lines(x=misclass90_Mantodea$Total,y=misclass90_Mantodea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[10])
lines(x=misclass90_Scolopendromorpha$Total,y=misclass90_Scolopendromorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[11])
lines(x=misclass90_Pantopoda$Total,y=misclass90_Pantopoda$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[12])
lines(x=misclass90_Amblypygi$Total,y=misclass90_Amblypygi$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[13])
lines(x=misclass90_Astigmata$Total,y=misclass90_Astigmata$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[14])
lines(x=misclass90_Blattodea$Total,y=misclass90_Blattodea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[15])
lines(x=misclass90_Cumacea$Total,y=misclass90_Cumacea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[16])
lines(x=misclass90_Dermaptera$Total,y=misclass90_Dermaptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[17])
lines(x=misclass90_Endeostigmata$Total,y=misclass90_Endeostigmata$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[18])
lines(x=misclass90_Glomerida$Total,y=misclass90_Glomerida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[19])
lines(x=misclass90_Grylloblattodea$Total,y=misclass90_Grylloblattodea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[20])
lines(x=misclass90_Halocryprida$Total,y=misclass90_Halocryprida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[21])
lines(x=misclass90_Harpacticoida$Total,y=misclass90_Harpacticoida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[22])
lines(x=misclass90_Julida$Total,y=misclass90_Julida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[23])
lines(x=misclass90_Lithobiomorpha$Total,y=misclass90_Lithobiomorpha$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[24])
lines(x=misclass90_Myodocopida$Total,y=misclass90_Myodocopida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[25])
lines(x=misclass90_Poecilostomatoida$Total,y=misclass90_Poecilostomatoida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[26])
lines(x=misclass90_Polyzoniida$Total,y=misclass90_Polyzoniida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[27])
lines(x=misclass90_Pseudoscorpiones$Total,y=misclass90_Pseudoscorpiones$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[28])
lines(x=misclass90_Psocoptera$Total,y=misclass90_Psocoptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[29])
lines(x=misclass90_Siphonophorida$Total,y=misclass90_Siphonophorida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[30])
lines(x=misclass90_Spirostreptida$Total,y=misclass90_Spirostreptida$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[31])
lines(x=misclass90_Strepsiptera$Total,y=misclass90_Strepsiptera$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[32])
lines(x=misclass90_Tanaidacea$Total,y=misclass90_Tanaidacea$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[33])
lines(x=misclass90_undef_Diplopoda$Total,y=misclass90_undef_Diplopoda$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[34])
lines(x=misclass90_undef_Malacostraca$Total,y=misclass90_undef_Malacostraca$Misclassified,type="p", pch=21, cex=1.25, col="black", bg=colors[35])

#x-axis
xticks<-c(1,10,100,1000,10000,100000)
axis(1, at=xticks,labels=xticks,pos=0)

#y-axis
#need two axis commands, one for the line, one for the ticks and labels
axis(2,at=c(0,100), labels=c("",""),lwd.ticks=0, pos=1)
axis(2, at=seq(0,100,by=20), lwd=0, lwd.ticks=1, pos=1, cex.axis=0.75)

lines(x=c(1,100000), y=c(10,10), lty=2, col="black")

par(mfrow=c(1,2))

#Create empty plot for legend only
plot(c(1,2,3), c(1,2,3), type="n", axes=FALSE, xlab="", ylab="")

legend(1, 3, legend=c("Raphidioptera, n=5", "Craterostigmomorpha, n=xx", "Siphonaptera, n=10", "Pauropoda, n=41", "Embioptera, n=54", "Hemiptera, n=8764", "Geophilomorpha, n=24", "Neuroptera, n=517", "Thysanoptera, n=696", "Mantodea, n=72", "Scolopendrmorpha, n=217", "Pantopoda, n=346", "Amblypygi, n=10", "Astigmata, n=78", "Blattodea, n=91", "Cumacea, n=24", "Dermaptera, n=30", "Endeostigmata, n=1", "Glomerida, n=8", "Grylloblattodea, n=1", "Halocryprida, n=44", "Harpacticoida, n=76", "Julida, n=27", "Lithobiomorpha, n=49", "Myodocopida, n=3", "Poecilostomatoida, n=41", "Polyzoniida, n=1", "Pseudoscorpiones, n=68", "Psocoptera, n=22", "Siphonophorida, n=2", "Spirostreptida, n=1", "Strepsiptera, n=22", "Tanaidacea, n=42", "undef_Diplopoda, n=9", "undef_Malacostraca, n=123"), cex=0.5, fill=colors, bty="n")

dev.off()

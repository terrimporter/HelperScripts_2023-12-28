library(vegan)

pdf("family_specaccum_plate.pdf")
par(mfrow=c(2,2),oma=c(0,0,2,0))
#title("Family", outer=TRUE)

#indvdenoised
A<-read.csv(file="family_abund.csv", head=TRUE, row.names=1)
A

#split out three datasets
A_NBE<-A[1:35,]
A_NBI<-A[36:71,]
A_NBR<-A[72:107,]

#remove columns with only zeros
A_NBE<-A_NBE[,colSums(A_NBI) !=0]
A_NBI<-A_NBI[,colSums(A_NBE) !=0]
A_NBR<-A_NBR[,colSums(A_NBR) !=0]

#run specaccum
A_NBE_specaccum<-specaccum(A_NBE)
A_NBI_specaccum<-specaccum(A_NBI)
A_NBR_specaccum<-specaccum(A_NBR)

#count sites
length_A_NBE_specaccum_sites<-length(A_NBE_specaccum$sites)
length_A_NBI_specaccum_sites<-length(A_NBI_specaccum$sites)
length_A_NBR_specaccum_sites<-length(A_NBR_specaccum$sites)

#globaldenoised
B<-read.csv(file="family_gOTU_abund.csv", head=TRUE, row.names=1)
B

#split out three datasets
B_NBE<-B[1:35,]
B_NBI<-B[36:71,]
B_NBR<-B[72:107,]

#remove columns with only zeros
B_NBE<-B_NBE[,colSums(B_NBI) !=0]
B_NBI<-B_NBI[,colSums(B_NBE) !=0]
B_NBR<-B_NBR[,colSums(B_NBR) !=0]

#run specaccum
B_NBE_specaccum<-specaccum(B_NBE)
B_NBI_specaccum<-specaccum(B_NBI)
B_NBR_specaccum<-specaccum(B_NBR)

#count sites
length_B_NBE_specaccum_sites<-length(B_NBE_specaccum$sites)
length_B_NBI_specaccum_sites<-length(B_NBI_specaccum$sites)
length_B_NBR_specaccum_sites<-length(B_NBR_specaccum$sites)

#indvdenoised97
C<-read.csv(file="family_abund_97.csv", head=TRUE, row.names=1)
C

#split out three datasets
C_NBE<-C[1:35,]
C_NBI<-C[36:71,]
C_NBR<-C[72:107,]

#remove columns with only zeros
C_NBE<-C_NBE[,colSums(C_NBI) !=0]
C_NBI<-C_NBI[,colSums(C_NBE) !=0]
C_NBR<-C_NBR[,colSums(C_NBR) !=0]

#run specaccum
C_NBE_specaccum<-specaccum(C_NBE)
C_NBI_specaccum<-specaccum(C_NBI)
C_NBR_specaccum<-specaccum(C_NBR)

#count sites
length_C_NBE_specaccum_sites<-length(C_NBE_specaccum$sites)
length_C_NBI_specaccum_sites<-length(C_NBI_specaccum$sites)
length_C_NBR_specaccum_sites<-length(C_NBR_specaccum$sites)

#globaldenoised97
D<-read.csv(file="family_gOTU_abund_97.csv", head=TRUE, row.names=1)
D

#split out three datasets
D_NBE<-D[1:35,]
D_NBI<-D[36:71,]
D_NBR<-D[72:107,]

#remove columns with only zeros
D_NBE<-D_NBE[,colSums(D_NBI) !=0]
D_NBI<-D_NBI[,colSums(D_NBE) !=0]
D_NBR<-D_NBR[,colSums(D_NBR) !=0]

#run specaccum
D_NBE_specaccum<-specaccum(D_NBE)
D_NBI_specaccum<-specaccum(D_NBI)
D_NBR_specaccum<-specaccum(D_NBR)

#count sites
length_D_NBE_specaccum_sites<-length(D_NBE_specaccum$sites)
length_D_NBI_specaccum_sites<-length(D_NBI_specaccum$sites)
length_D_NBR_specaccum_sites<-length(D_NBR_specaccum$sites)

#find max richness across all 4 plots
combinedrichness<-c(A_NBE_specaccum$richness, A_NBI_specaccum$richness, A_NBR_specaccum$richness, B_NBE_specaccum$richness, B_NBI_specaccum$richness, B_NBR_specaccum$richness,C_NBE_specaccum$richness, C_NBI_specaccum$richness, C_NBR_specaccum$richness, D_NBE_specaccum$richness, D_NBI_specaccum$richness, D_NBR_specaccum$richness)
maxrichness<-max(combinedrichness)

#find max sites across all 4 plots
combinedsites<-c(length_A_NBE_specaccum_sites, length_A_NBI_specaccum_sites, length_A_NBR_specaccum_sites, length_B_NBE_specaccum_sites, length_B_NBI_specaccum_sites, length_B_NBR_specaccum_sites, length_C_NBE_specaccum_sites,length_C_NBI_specaccum_sites, length_C_NBR_specaccum_sites, length_D_NBE_specaccum_sites,length_D_NBI_specaccum_sites, length_D_NBR_specaccum_sites)
combinedsites
maxsites<-max(combinedsites)
maxsites

#three curves in one plot
plot(range(c(1,maxsites)), range(c(1,maxrichness)), xlab="Samples", ylab="No. families", main="Indv Denoised", type="n")

lines(A_NBE_specaccum$sites, A_NBE_specaccum$richness, lwd=1.5, lty=1, col="red4")
upper<-mapply("+",A_NBE_specaccum$richness, A_NBE_specaccum$sd)
lines(A_NBE_specaccum$sites, upper, lwd=0.5, lty=2, col="red4")
lower<-mapply("-",A_NBE_specaccum$richness, A_NBE_specaccum$sd)
lines(A_NBE_specaccum$sites, lower, lwd=0.5, lty=2, col="red4")

lines(A_NBI_specaccum$sites, A_NBI_specaccum$richness, lwd=1.5, lty=1, col="blue4")
upper<-mapply("+",A_NBI_specaccum$richness, A_NBI_specaccum$sd)
lines(A_NBI_specaccum$sites, upper, lwd=0.5, lty=2, col="blue4")
lower<-mapply("-",A_NBI_specaccum$richness, A_NBI_specaccum$sd)
lines(A_NBI_specaccum$sites, lower, lwd=0.5, lty=2, col="blue4")

lines(A_NBR_specaccum$sites, A_NBR_specaccum$richness, lwd=1.5, lty=1, col="green4")
upper<-mapply("+",A_NBR_specaccum$richness, A_NBR_specaccum$sd)
lines(A_NBR_specaccum$sites, upper, lwd=0.5, lty=2, col="green4")
lower<-mapply("-",A_NBR_specaccum$richness, A_NBR_specaccum$sd)
lines(A_NBR_specaccum$sites, lower, lwd=0.5, lty=2, col="green4")

legend("bottomright", legend=c("NBI","NBE","NBR"), col=c("red4","blue4","green4"), lty=1, lwd=1.5, bty="n")

#three curves in one plot
plot(range(c(1,maxsites)), range(c(1,maxrichness)), xlab="Samples", ylab="No. families", main="Global Denoised", type="n")

lines(B_NBE_specaccum$sites, B_NBE_specaccum$richness, lwd=1.5, lty=1, col="red4")
upper<-mapply("+",B_NBE_specaccum$richness, B_NBE_specaccum$sd)
lines(B_NBE_specaccum$sites, upper, lwd=0.5, lty=2, col="red4")
lower<-mapply("-",B_NBE_specaccum$richness, B_NBE_specaccum$sd)
lines(B_NBE_specaccum$sites, lower, lwd=0.5, lty=2, col="red4")

lines(B_NBI_specaccum$sites, B_NBI_specaccum$richness, lwd=1.5, lty=1, col="blue4")
upper<-mapply("+",B_NBI_specaccum$richness, B_NBI_specaccum$sd)
lines(B_NBI_specaccum$sites, upper, lwd=0.5, lty=2, col="blue4")
lower<-mapply("-",B_NBI_specaccum$richness, B_NBI_specaccum$sd)
lines(B_NBI_specaccum$sites, lower, lwd=0.5, lty=2, col="blue4")

lines(B_NBR_specaccum$sites, B_NBR_specaccum$richness, lwd=1.5, lty=1, col="green4")
upper<-mapply("+",B_NBR_specaccum$richness, B_NBR_specaccum$sd)
lines(B_NBR_specaccum$sites, upper, lwd=0.5, lty=2, col="green4")
lower<-mapply("-",B_NBR_specaccum$richness, B_NBR_specaccum$sd)
lines(B_NBR_specaccum$sites, lower, lwd=0.5, lty=2, col="green4")

#three curves in one plot
plot(range(c(1,maxsites)), range(c(1,maxrichness)), xlab="Samples", ylab="No. families", main="Indv Denoised 97%", type="n")

lines(C_NBE_specaccum$sites, C_NBE_specaccum$richness, lwd=1.5, lty=1, col="red4")
upper<-mapply("+",C_NBE_specaccum$richness, C_NBE_specaccum$sd)
lines(C_NBE_specaccum$sites, upper, lwd=0.5, lty=2, col="red4")
lower<-mapply("-",C_NBE_specaccum$richness, C_NBE_specaccum$sd)
lines(C_NBE_specaccum$sites, lower, lwd=0.5, lty=2, col="red4")

lines(C_NBI_specaccum$sites, C_NBI_specaccum$richness, lwd=1.5, lty=1, col="blue4")
upper<-mapply("+",C_NBI_specaccum$richness, C_NBI_specaccum$sd)
lines(C_NBI_specaccum$sites, upper, lwd=0.5, lty=2, col="blue4")
lower<-mapply("-",C_NBI_specaccum$richness, C_NBI_specaccum$sd)
lines(C_NBI_specaccum$sites, lower, lwd=0.5, lty=2, col="blue4")

lines(C_NBR_specaccum$sites, C_NBR_specaccum$richness, lwd=1.5, lty=1, col="green4")
upper<-mapply("+",C_NBR_specaccum$richness, C_NBR_specaccum$sd)
lines(C_NBR_specaccum$sites, upper, lwd=0.5, lty=2, col="green4")
lower<-mapply("-",C_NBR_specaccum$richness, C_NBR_specaccum$sd)
lines(C_NBR_specaccum$sites, lower, lwd=0.5, lty=2, col="green4")

#three curves in one plot
plot(range(c(1,maxsites)), range(c(1,maxrichness)), xlab="Samples", ylab="No. families", main="Global Denoised 97%", type="n")

lines(D_NBE_specaccum$sites, D_NBE_specaccum$richness, lwd=1.5, lty=1, col="red4")
upper<-mapply("+",D_NBE_specaccum$richness, D_NBE_specaccum$sd)
lines(D_NBE_specaccum$sites, upper, lwd=0.5, lty=2, col="red4")
lower<-mapply("-",D_NBE_specaccum$richness, D_NBE_specaccum$sd)
lines(D_NBE_specaccum$sites, lower, lwd=0.5, lty=2, col="red4")

lines(D_NBI_specaccum$sites, D_NBI_specaccum$richness, lwd=1.5, lty=1, col="blue4")
upper<-mapply("+",D_NBI_specaccum$richness, D_NBI_specaccum$sd)
lines(D_NBI_specaccum$sites, upper, lwd=0.5, lty=2, col="blue4")
lower<-mapply("-",D_NBI_specaccum$richness, D_NBI_specaccum$sd)
lines(D_NBI_specaccum$sites, lower, lwd=0.5, lty=2, col="blue4")

lines(D_NBR_specaccum$sites, D_NBR_specaccum$richness, lwd=1.5, lty=1, col="green4")
upper<-mapply("+",D_NBR_specaccum$richness, D_NBR_specaccum$sd)
lines(D_NBR_specaccum$sites, upper, lwd=0.5, lty=2, col="green4")
lower<-mapply("-",D_NBR_specaccum$richness, D_NBR_specaccum$sd)
lines(D_NBR_specaccum$sites, lower, lwd=0.5, lty=2, col="green4")

title("Family", outer=TRUE)

dev.off()

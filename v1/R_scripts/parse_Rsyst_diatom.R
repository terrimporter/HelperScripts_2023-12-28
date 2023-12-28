# Teresita M. Porter, Feb. 11/20
# Script to parse R-syst diatom release v7.1 into a fasta file that can be parsed into files for the RDP classifier

A <- read.csv("2019-02-12-Diat.barcode-release-version 7.1.csv", header = TRUE, stringsAsFactors = FALSE)

# subset by amplicon
SSU <- A[A$Amplified.region=="18S" |
         A$Amplified.region=="18S-V4" |
         A$Amplified.region=="18s" |
         A$Amplified.region=="18s-V4" |
         A$Amplified.region=="28S",]
# 2962

rbcL <- A[A$Amplified.region=="rbcl" |
          A$Amplified.region=="Rbcl-312bp",]
# 3504

# subset to keep id, lineage, seq
SSU.1 <- SSU[,c(1,3,36,35,34,33,32,31,30,29,28)]

rbcL.1 <- rbcL[,c(1,3,36,35,34,33,32,31,30,29,28)]

# save to csv file
write.table(SSU.1, file="SSU.csv", sep=",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(rbcL.1, file="rbcL.csv", sep=",", col.names = TRUE, row.names = FALSE, quote = FALSE)

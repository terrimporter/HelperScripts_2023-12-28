library(vegan)
library(ggplot2)
library(ggpubr)
library(reshape)
library(splitstackshape)
library(scales)

# Read infile prepared by python script
A<-read.table(file="ESV_LabCode_SiteCode_matrix.csv", head=TRUE, row.names=1, sep=",")
str(A)

# Create box plot of reads vs labs & sites

## Sum total number of reads per sample
LabCodeReads<-rowSums(A)

## Rename columns
df<-data.frame("Reads"=LabCodeReads)
df$LabCode_SiteCode<-rownames(df)

## Copy a column for splitting
df$combo<-df$LabCode_SiteCode

## Split combo into new ones
df2 <- cSplit(df,"combo", sep="_", direction="wide", fixed=T)

## Rename new columns
names(df2)[3]<-"LabCode"
names(df2)[4]<-"SiteCode"

## Copy column
df2$Site<-df2$SiteCode

## Rename numbered sites into their site names
df2$Site <- gsub('10', 'Hillary', df2$Site)
df2$Site <- gsub('11', 'Hillary', df2$Site)
df2$Site <- gsub('12', 'Hillary', df2$Site)
df2$Site <- gsub('13', 'Waitemata', df2$Site)
df2$Site <- gsub('14', 'Waitemata', df2$Site)
df2$Site <- gsub('15', 'Waitemata', df2$Site)
df2$Site <- gsub('16', 'Waitemata', df2$Site)
df2$Site <- gsub('17', 'Blank', df2$Site)
df2$Site <- gsub('18', 'Blank', df2$Site)
df2$Site <- gsub('19', 'NoTemplatePositive', df2$Site)
df2$Site <- gsub('20', 'NoTemplatePositive', df2$Site)
df2$Site <- gsub('1', 'Monterey', df2$Site)
df2$Site <- gsub('2', 'Monterey', df2$Site)
df2$Site <- gsub('3', 'Monterey', df2$Site)
df2$Site <- gsub('4', 'Monterey', df2$Site)
df2$Site <- gsub('5', 'Victoria', df2$Site)
df2$Site <- gsub('6', 'Victoria', df2$Site)
df2$Site <- gsub('7', 'Victoria', df2$Site)
df2$Site <- gsub('8', 'Victoria', df2$Site)
df2$Site <- gsub('9', 'Hillary', df2$Site)

## Exclude Blank and NoTemplatePositive rows
df3 <- subset(df2, Site!="Blank")
df4 <- subset(df3, Site!="NoTemplatePositive")

## Box plot reads vs lab codes, colored by site
Reads_boxplot<-ggplot(df4, aes(x = LabCode, y = Reads, fill = Site)) + 
  ggtitle("Uni18S + Uni18SR") +
  ylab("Denoised reads") +
  geom_boxplot() +
  scale_y_continuous(labels=comma)
ggsave("Reads.pdf")

# Create box plot for ESVs vs lab & site

## Transform from abundance to presence/absence
A[A > 0] <- 1 

## Sum total number of reads per sample
LabCodeESVs<-rowSums(A)

## Rename columns
df5<-data.frame("ESVs"=LabCodeESVs)
df5$LabCode_SiteCode<-rownames(df5)

## Copy a column for splitting
df5$combo<-df5$LabCode_SiteCode

## Split combo into new ones
df6 <- cSplit(df5,"combo", sep="_", direction="wide", fixed=T)

## Rename new columns
names(df6)[3]<-"LabCode"
names(df6)[4]<-"SiteCode"

## Copy column
df6$Site<-df6$SiteCode

## Rename numbered sites into their site names
df6$Site <- gsub('10', 'Hillary', df6$Site)
df6$Site <- gsub('11', 'Hillary', df6$Site)
df6$Site <- gsub('12', 'Hillary', df6$Site)
df6$Site <- gsub('13', 'Waitemata', df6$Site)
df6$Site <- gsub('14', 'Waitemata', df6$Site)
df6$Site <- gsub('15', 'Waitemata', df6$Site)
df6$Site <- gsub('16', 'Waitemata', df6$Site)
df6$Site <- gsub('17', 'Blank', df6$Site)
df6$Site <- gsub('18', 'Blank', df6$Site)
df6$Site <- gsub('19', 'NoTemplatePositive', df6$Site)
df6$Site <- gsub('20', 'NoTemplatePositive', df6$Site)
df6$Site <- gsub('1', 'Monterey', df6$Site)
df6$Site <- gsub('2', 'Monterey', df6$Site)
df6$Site <- gsub('3', 'Monterey', df6$Site)
df6$Site <- gsub('4', 'Monterey', df6$Site)
df6$Site <- gsub('5', 'Victoria', df6$Site)
df6$Site <- gsub('6', 'Victoria', df6$Site)
df6$Site <- gsub('7', 'Victoria', df6$Site)
df6$Site <- gsub('8', 'Victoria', df6$Site)
df6$Site <- gsub('9', 'Hillary', df6$Site)

## Exclude Blank and NoTemplatePositive rows
df7 <- subset(df6, Site!="Blank")
df8 <- subset(df7, Site!="NoTemplatePositive")

## Create box plot reads vs lab codes, colored by site
ESVs_boxplot<-ggplot(df8, aes(x = LabCode, y = ESVs, fill = Site)) + 
  ggtitle("Uni18S + Uni18SR") +
  ylab("ESVs") +
  geom_boxplot() +
  scale_y_continuous(labels=comma)
ggsave("ESVs.pdf")

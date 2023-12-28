library(vegan)
library(ggplot2)
library(ggpubr)
library(reshape)
library(splitstackshape)
library(scales)

# Create stacked bar plots for target taxa, reads vs lab group, sites pooled

## Read infile prepared by python script
A<-read.table(file="AllPhyla_LabCode_SiteCode_matrix.csv", head=TRUE, row.names=1, sep=",")
#str(A)
#head(A)

## Transpose
At<-t(A)

## Reshape for ggplot
long<-melt(At)
head(long)

## Rename new columns
names(long)[1]<-"Phyla"
names(long)[2]<-"LabCode_SiteCode"
names(long)[3]<-"Reads"

## Split LabCode_SiteCode into their own columns
long2 <- cSplit(long,"LabCode_SiteCode", sep="_", direction="wide", fixed=T)

## Rename new columns
names(long2)[1]<-"Phyla"
names(long2)[2]<-"Reads"
names(long2)[3]<-"LabCode"
names(long2)[4]<-"SiteCode"
head(long2)

## Copy column then add actual site names
long2$Site<-long2$SiteCode

## Rename numbered sites into their site names
long2$Site <- gsub('10', 'Hillary', long2$Site)
long2$Site <- gsub('11', 'Hillary', long2$Site)
long2$Site <- gsub('12', 'Hillary', long2$Site)
long2$Site <- gsub('13', 'Waitemata', long2$Site)
long2$Site <- gsub('14', 'Waitemata', long2$Site)
long2$Site <- gsub('15', 'Waitemata', long2$Site)
long2$Site <- gsub('16', 'Waitemata', long2$Site)
long2$Site <- gsub('17', 'Blank', long2$Site)
long2$Site <- gsub('18', 'Blank', long2$Site)
long2$Site <- gsub('19', 'NoTemplatePositive', long2$Site)
long2$Site <- gsub('20', 'NoTemplatePositive', long2$Site)
long2$Site <- gsub('1', 'Monterey', long2$Site)
long2$Site <- gsub('2', 'Monterey', long2$Site)
long2$Site <- gsub('3', 'Monterey', long2$Site)
long2$Site <- gsub('4', 'Monterey', long2$Site)
long2$Site <- gsub('5', 'Victoria', long2$Site)
long2$Site <- gsub('6', 'Victoria', long2$Site)
long2$Site <- gsub('7', 'Victoria', long2$Site)
long2$Site <- gsub('8', 'Victoria', long2$Site)
long2$Site <- gsub('9', 'Hillary', long2$Site)

## Exclude Blank and NoTemplatePositive rows
long3 <- subset(long2, Site!="Blank")
long4 <- subset(long3, Site!="NoTemplatePositive")

## Create stacked bar plot reads vs lab codes, colored by phyla
Phyla_stacked<-ggplot(long4) + 
  ggtitle("Uni18S + Uni18SR") +
  ylab("Reads") +
  geom_bar(aes(x=LabCode, y=Reads,fill=Phyla), color="black", stat="identity") +
  scale_y_continuous(labels=comma) +
  guides(fill=guide_legend(title="All Eukaryota Phyla",
			                     keywidth=0.1,
                            keyheight=0.1,
                          default.unit="inch",
			                     title.position="top"),
		          shape = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 7),
		  legend.position = "bottom")
ggsave("StackedBar_AllPhyla.pdf")

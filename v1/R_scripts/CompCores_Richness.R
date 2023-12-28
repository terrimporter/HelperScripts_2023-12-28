#ILC_richness.csv and NZC_richness.csv compes from CompareCores_ESVs_genera_families_richness.py

###################################################################
##### Richness bar plot
###################################################################

#read infile prepared by python script
#contains richness summarized for ESVs, genera, families
#1C samples were subsampeld to N=4 to match S1-S4 for XC samples
I_richness<-read.csv(file="ILC_richness.csv", head=TRUE, row.names=1)
I_richness

N_richness<-read.csv(file="NZC_richness.csv", head=TRUE, row.names=1)
N_richness

library(ggplot2)
library(reshape)

#reshape ILC data for ggplot
I_richness2<-melt(I_richness, 
                id.vars=c("rank"))
I_richness2

Ig <- ggplot(I_richness2, aes(fill=variable,y=value,x=rank)) +
  geom_bar(position="dodge",stat="identity") +
  theme_bw() +
  ggtitle("ILC45") +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        plot.title=element_text(size=20, hjust=0.5)) +
  xlab("Rank") +
  ylab("Richness (log10)") +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=3) +
  scale_fill_discrete(name="Cores",
                      breaks=c("X1C","X2C","X4C","X6C","X8C","12_15C"),
                      labels=c("1","2","4","6","8","12-15"))
  
pdf("ILC_richness.pdf")
Ig
dev.off()

#reshape NZC data for ggplot
N_richness2<-melt(N_richness, 
                  id.vars=c("rank"))
N_richness2

Ng <- ggplot(N_richness2, aes(fill=variable,y=value,x=rank)) +
  geom_bar(position="dodge",stat="identity") +
  theme_bw() +
  ggtitle("NZC85") +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        plot.title=element_text(size=20, hjust=0.5)) +
  xlab("Rank") +
  ylab("Richness (log10)") +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=3) +
  scale_fill_discrete(name="Cores",
                      breaks=c("X1C","X2C","X4C","X6C","X8C","12_15C"),
                      labels=c("1","2","4","6","8","12-15"))

pdf("NZC_richness.pdf")
Ng
dev.off()

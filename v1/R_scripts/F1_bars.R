library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)

# CSV file from excel
# Figure 1a) year vs Deposited
# Figure 1b) year vs UniqueSpecies
# Figure 1c) year vs Deposited (log) for 2003 to 2003-2017 only

# Read in file
df<-read.csv("F1.csv",header=TRUE)

# Remove 2003_2017 (last row) 
df2<-df[-nrow(df),]

# Deposited
df3<-df2[,1:2]

# UniqueSpecies
df4<-df2[,c(1,3)]

# 2003 vs 2003-2017
df5<-df[c(15,16),]

# Melt so two series can be plotted, id.vars=Year, fix legend labels here
df5.long<-melt(df5, id.vars="Year", variable.name="Series", value.name="Count")
df5.long$Series<-factor(df5.long$Series, 
                        levels=c("Deposited","UniqueSpecies"), 
                        labels=c("Deposited","Unique Species"))

# create two bar plots
F1a<-ggplot(df3, aes(Year)) +
  ggtitle("a)\n") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        axis.text.x = element_text(angle=90)) +
  geom_bar(aes(weight=Deposited), color="#F8766D", fill="#F8766D", width=0.8) +
  xlab("Year") +
  ylab("Records") +
  scale_y_continuous(label=comma)

F1b<-ggplot(df4, aes(Year)) +
  ggtitle("b)\n") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        axis.text.x = element_text(angle=90)) +
  geom_bar(aes(weight=UniqueSpecies), color="#00BFC4", fill="#00BFC4", width=0.8) +
  xlab("Year") +
  ylab("Unique Species") +
  scale_y_continuous(label=comma)

# Create line plot with two series
F1c<-ggplot(df5.long,aes(x=Year,y=Count, color=Series, group=Series, shape=Series)) +
  theme_bw() +
  ggtitle("c)\n") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        #margin order start at top, clockwise
        plot.margin=margin(0,0,0,1,"cm")) +
  guides(colour=guide_legend(nrow=2)) +
  geom_point(size=3) +
  geom_text(aes(label=comma(Count)),hjust=0,vjust=1, nudge_x=0.1, nudge_y=0.1, color="black", size=4) +
  geom_line() +
  xlab("Year") +
  ylab("Number (Log10)") +
  #expand adds spacing around plotted points so the text labels aren't cutoff
  scale_x_discrete(labels=c("2003"="2003","2003_2017"="2003-2017"), expand = c(1,1)) +
  scale_y_log10(label=comma)

#define subplot layout
lay<-rbind(c(1,3),
           c(2,3))

pdf("Fig1.pdf")
grid.arrange(F1a,F1b,F1c, layout_matrix=lay)
dev.off()
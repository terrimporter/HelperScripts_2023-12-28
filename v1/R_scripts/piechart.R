orders <- read.table(file="iBOL_all.csv",header=FALSE,sep=",")
orders
is.list(orders)

slices <- orders$V2
slices

lbls <- orders$V1
lbls

pct <- round(slices/sum(slices)*100)

lbls <- paste(lbls, slices, sep="\n")
lbls <- paste(lbls, pct, sep="\n")
lbls <- paste(lbls, "%", sep="")

library(RColorBrewer)
#need 11 colors
colors1 = brewer.pal(9,"Set1")
colors2 = brewer.pal(2,"Set2")
colors = c(colors1, colors2)

pie(slices, labels=lbls, main="iBOL data release 3.75 v1\nTaxonomic order abundance", col=colors, clockwise=TRUE)

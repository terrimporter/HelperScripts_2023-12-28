# Teresita M. Porter, Dec. 13, 2019

df <- read.table("lengths.txt", header=FALSE)

lq <- quantile(df[,1], 0.25) 
uq <- quantile(df[,1], 0.75)

iqr <- uq - lq

min <- lq - (1.5*iqr)
max <- uq + (1.5*iqr)

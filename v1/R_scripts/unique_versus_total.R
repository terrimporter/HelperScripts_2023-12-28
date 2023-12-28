a<-read.table("incidence.csv",sep=",",header=TRUE)
#remove sites column
#a$X
b<-a[,-1]
#colsums<-colSums(b)
rowsums<-rowSums(b)
rowsums

#get number of sites in rows
#length_b<-length(b[,1])

#row<-vector()
#counter<-0
#unique<-0
#percents<-vector()

#for(i in c(1:length_b)){
#	row<-b[i,]
#	for(j in row){
#		counter=counter+1
#		if(row[counter]==colsums[counter]){
#			unique=unique+1
#		}
#	}
#	counter<-0
#	sum<-rowSums(b[i,])
#	percent<-(unique/sum)*100
#	percents<-append(percents,percent)
#	unique<-0
#	row<-vector()
#}
#percents

a<-read.table("incidence.csv",sep=",",header=TRUE)
#remove sites column
b<-a[,-1]

#total number of OTUs in subsampled incidence.csv
number_columns<-ncol(b)
number_columns

#total per OTU
colsums<-colSums(b)
colsums

#total per subsite
rowsums<-rowSums(b)
rowsums

#get number of subsites in rows
length_b<-length(b[,1])

row<-vector()
counter<-0
unique<-0
uniques<-vector()
sums<-vector()
percents<-vector()

#for each row (subsite)
for(i in 1:length_b){
	row<-b[i,]

	for(j in row){
		counter=counter+1
#		print(row[counter])
#		print(colsums[counter])
		if(row[counter]>0) {
			if(row[counter]==colsums[counter]){
				unique=unique+1
			}
		}
	}
	uniques<-append(uniques,unique)
	sum<-rowsums[i]
	sums<-append(sums,sum)
	percent<-(unique/sum)*100
	percents<-append(percents,percent)
	
	row<-vector()
	counter<-0
	unique<-0
}

percents

#create df for just PAD
length_b<-length(b[1,]) #get total number of columns

pad<-b[1:8,]
#pad
colsums<-colSums(pad)
colsums 

length_pad<-length(pad[1,]) #number of columns

counter<-0

for (i in 1:length_pad){
	if (colsums[i] > 0) {
#		PAD[i]
		counter = counter+1
	}
}
#number of OTUs in PAD
counter

#create df for just WC
WC<-b[9:16,]
colsums<-colSums(WC)
colsums

counter<-0

for (i in 1:length(colsums)) {
	if (colsums[i] > 0) {
		counter = counter+1
	}
}
#number of OTUs in WC
counter

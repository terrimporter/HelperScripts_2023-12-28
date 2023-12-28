library(vegan)
pdf("1C1E_rarecurve.pdf")
par(mfrow=c(1,1), mar=c(7,4,4,2))

#indvdenoised_97OTUs
A<-read.csv(file="1C1E_OTU_matrix.csv", head=TRUE, row.names=1)
#A

#split out three datasets
#A_part<-A[1:48,]

#remove columns with only zeros
A_notnull<-A[,colSums(A) !=0]
#A_notnull<-A_part[,colSums(is.na(A_part)) != nrow(A_part)]
#A_notnull

#remove rows with only zeros
A_notnull2<-A_notnull[rowSums(A_notnull) !=0,]

#get total number of taxa per site
S<-specnumber(A_notnull2)
S
#get largest library size with most reads
readmax<-max(rowSums(A_notnull2))
readmax
#get smallest library size with least reads
readmin<-min(rowSums(A_notnull2))
readmin
#get largest library richness with most taxa
taxmax<-max(S)
taxmax
#edit vegan function rarecurve to fix plot labels
rarecurve_ESV<-function (x, step, sample,main="2016 Grid 1C1E Expt",xlab="Reads", ylab="ESVs",label=TRUE, col, lty, ...){
	x <- as.matrix(x)
	if (!identical(all.equal(x, round(x)), TRUE))
	stop("function accepts only integers (counts)")
	if (missing(col))
	col <- par("col")
	if (missing(lty))
	lty <- par("lty")
	tot <- rowSums(x)
	S <- specnumber(x)
	if (any(S <= 0)) {
		message("empty rows removed")
		x <- x[S > 0, , drop = FALSE]
		tot <- tot[S > 0]
		S <- S[S > 0]
	}
	nr <- nrow(x)
	col <- rep(col, length.out = nr)
	lty <- rep(lty, length.out = nr)
	out <- lapply(seq_len(nr), function(i) {
		n <- seq(1, tot[i], by = step)
		if (n[length(n)] != tot[i])
		n <- c(n, tot[i])
		drop(rarefy(x[i, ], n))
	})
	Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
	Smax <- sapply(out, max)
	plot(c(1, max(Nmax)), c(1, max(Smax)),main=main, xlab = xlab, ylab = ylab,type = "n", ...)
	if (!missing(sample)) {
		abline(v = sample)
		rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),y = z, xout = sample, rule = 1)$y)
		abline(h = rare, lwd = 0.5)
	}
	for (ln in seq_along(out)) {
		N <- attr(out[[ln]], "Subsample")
		lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)						    			   }
	if (label) {
		ordilabel(cbind(tot, S), labels = rownames(x), ...)
	}
	invisible(out)
}

#add third var 'sample=readmin' to see horizontal lines that show rarified richness for each sampleat this depth of sampling.  Ignored here because smallest library size is just 5 (after removing singletons and doubletons)
rarecurveout<-rarecurve_ESV(A_notnull2, step=100, col="blue4", cex=0.6)
readmin
dev.off()

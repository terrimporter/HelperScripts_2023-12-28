library(vegan)
pdf("family_rarecurve.pdf")
par(mfrow=c(1,1), mar=c(7,4,4,2))

#indvdenoised_97OTUs
A<-read.csv(file="family_abund_denoised_97.csv", head=TRUE, row.names=1)
#A

#split out three datasets
A_DT<-A[1:107,]

#remove columns with only zeros
A_DT<-A_DT[,colSums(A_DT) !=0]

#get total number of taxa per site
S<-specnumber(A_DT)

#get largest library size with most reads
readmax<-max(rowSums(A_DT))

#get smallest library size with least reads
readmin<-min(rowSums(A_DT))

#get largest library richness with most taxa
taxmax<-max(S)

#edit vegan function rarecurve to fix plot labels
#red=Tube03
#blue=Tube04
rarecurve_family<-function (x, step, sample,main="Greer01_Tube03_Tube04 DK ITS",xlab = "Reads", ylab = "Families",label = TRUE, 
	cols=c('red4','red4','red4','red4','blue4',
		'red4','red4','red4','red4','red4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',

		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',

		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4',
		'red4','red4','red4','red4','blue4','blue4'),
	lty, ...){
	x <- as.matrix(x)
	if (!identical(all.equal(x, round(x)), TRUE))
	stop("function accepts only integers (counts)")
#	if (missing(col))
#	col <- par("col")
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
#edit to plot lines in different colors
#	col <- rep(col, length.out = nr)
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
	legend("bottomright", legend=c("Tube03","Tube04"), col=c("red4","blue4"), lty=1, lwd=1.5, bty="n")
	if (!missing(sample)) {
		abline(v = sample)
		rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),y = z, xout = sample, rule = 1)$y)
		abline(h = rare, lwd = 0.5)
	}
	for (ln in seq_along(out)) {
		N <- attr(out[[ln]], "Subsample")
		#edit to show plotted curves in different colors
		lines(N, out[[ln]], col = cols[ln], lty = lty[ln], ...)						    			   }
	if (label) {
		ordilabel(cbind(tot, S), labels = rownames(x), ...)
	}
	invisible(out)
}

#add third var 'sample=readmin' to see horizontal lines that show rarified richness for each sampleat this depth of sampling.  Ignored here because smallest library size is just 5 (after removing singletons and doubletons)
rarecurveout<-rarecurve_family(A_DT, step=100, cex=0.6)
readmin
dev.off()


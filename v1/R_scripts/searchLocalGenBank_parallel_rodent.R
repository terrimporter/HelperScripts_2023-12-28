# Teresita M. Porter, Dec. 3, 2019

library(restez) # to search local GenBank database
#install.packages("data.table")
require(data.table) # fread function to read each filename in a list
library(readr) # read string from file
library(stringr) # regex capture group
#install.packages("read.gb")
#install.packages("seqinr")
library(seqinr) # complement
#install.packages("foreach")
library(foreach)
#install.packages("doParallel")
library(doParallel)

################################################################################
# function to get substring, complement, add to list
################################################################################
substring_comp <- function (gene.location,record) {
  result <- str_match(gene.location, "(\\d+)>?\\.\\.(\\d+)")
  start <- result[,2]
  stop <- result[,3]
  seq.all <- gb_extract(record = record, what = 'sequence')
  seq.sub <- substr(seq.all, start, stop)
  seq.sub.v <- s2c(seq.sub)
  # complement seq
  seq.v <- comp(seq.sub.v, forceToLower = FALSE, ambiguous = TRUE)
  seq <- paste(seq.v, collapse='')
  return.results<-c(seq,start,stop)
  return(return.results)
}
################################################################################
# function to parse a record
################################################################################
parse_record <- function (record, target, i, gb.seq, gb.taxid) {
  
#  record <- recs[1] #test
  accession <- gb_extract(record = record, what = 'accession')
  ft <- gb_extract(record = record, what = 'features')
  
  for (j in 1:length(ft)) {
    # parse through source part of features
    if (is.character(ft[[j]]$db_xref)) {
      taxid.field <- ft[[j]]$db_xref
      taxon.field <- strsplit(taxid.field, "\\s")[[1]][1] # sometimes coord here too <1..>234
      taxid <- str_match(taxon.field, "(\\d+)")[,1]
    }
    if (is.character(ft[[j]]$gene)) {
      gene.field <- ft[[j]]$gene
      if (gene.field %in% target) {
        # keep whole COI including introns (differs from previous versions of classifier <= v4)
        gene.location <- ft[[j]]$location
        if (grepl("complement",gene.location)) {
          return.results <- substring_comp(gene.location, record) # custom function above
          seq<-return.results[1]
          start<-return.results[2]
          stop<-return.results[3]
          gb.seq[[i]] <- paste(accession,seq, sep=" ")
          gb.taxid[[i]] <- paste(accession,taxid, sep=" ")
        }
        else {
          result <- str_match(gene.location, "<?(\\d+)>?\\.\\.(\\d+)")
          start <- result[,2]
          stop <- result[,3]
          seq.all <- gb_extract(record = record, what = 'sequence')[[1]]
          seq <- substr(seq.all, start, stop)
          gb.seq[[i]] <- paste(accession,seq, sep=" ")
          gb.taxid[[i]] <- paste(accession,taxid, sep=" ")
        }
      }
    }
  }
  map.list <- list(gb.seq, gb.taxid)
  return(map.list)
}
################################################################################
# function that transposes the results
# https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
################################################################################

# set the restez path to a memorable location
rodents_path = "/home/terri/rodents"
restez_path_set(rodents_path)

# download for domain 16 'rodents' # done
db_download(preselection = '16')

# create the rodent dbase for testing only
db_create(min_length = 100, max_length = 1000)

# Connect to SQL database in restez path
restez_connect()

# check status
restez_status()

# Search NCBI nucleotide database, retrieve accessions (small, fast)

# define gene names to target cytochrome oxidase c, subunit I in Eukaryotes
genes <- '("cox1"[GENE] OR "coxI"[GENE] OR "CO1"[GENE] OR "COI"[GENE])'
# minimum length 500 bp
# 99999999 indicates infinity
length <- "500:99999999[SLEN]"
# up to and including 2019
# use 2 as the lower bound
#date <- "2:2019[PDAT]"

# initialize list
accessions_cat <- list()

# read string of taxa from file
search_term <- paste(genes, "AND", length, "AND Rodentia[ORGN]", sep=" ")

search_object <- rentrez::entrez_search(db = 'nucleotide', 
	                                      term = search_term, 
	                                      use_history=TRUE, 
	                                      retmax = 0)

# searches local database first, then GenBank online
accession.line <- rentrez::entrez_fetch(db = 'nucleotide',
    	                        web_history = search_object$web_history,
        	                    rettype = 'acc')

# separate onto their own lines
accessions_cat <- strsplit(x = accession.line, split = '\n')[[1]]

# get the GenBank records from local rodent database
# a vector of records
recs <- gb_record_get(id=accessions_cat)

# specify target gene names
target <- c("cox1","coxI","co1","coI","COX1","COXI","CO1","COI")

# initialize lists
gb.taxid <- list()
gb.seq <- list()

# setup parallel backend to use many processors
# https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r
cores=12
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# parse records in parallel
finalList <- foreach(i=1:length(recs), .combine="comb", .multicombine=TRUE, .init=list(list(),list()), .packages = c("readr", "restez", "stringr", "seqinr")) %dopar% {
  record <- recs[i]
  map.list <- parse_record(record, target, i, gb.seq, gb.taxid)
  gb.seq <- map.list[[1]] # gb.seq
  gb.taxid <- map.list[[2]] # gb.taxid
  list(gb.seq, gb.taxid)
}

#stop cluster
stopCluster(cl)

# extract each list, turn into vector, recursive=TRUE (default)
gb.seq <- unlist(finalList[[1]])
gb.taxid <- unlist(finalList[[2]])

# print out the two mapping files
write.table(gb.taxid,"/home/terri/CO1v5/test/gb_taxid.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(gb.seq,"/home/terri/CO1v5/test/gb_seq.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Disconnect
restez_disconnect()







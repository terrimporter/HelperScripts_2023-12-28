# install.packages("devtools")
#library(devtools)
# install key dependency no longer available on CRAN
#Sys.setenv(TAR = "/bin/tar")
#devtools::install_github("hannesmuehleisen/MonetDBLite-R")
# install restez
#devtools::install_github("ropensci/restez")
library(restez)


# set path to local GenBank database
db <- "/home/terri/GenBank"
restez_path_set(filepath=db)

# download GenBank files (skip Bacteria for now) SLOW BE PATIENT!
db_download(preselection = '6 7 9 12 14 16')
# 6 Plants, fungi, algae
# 7 Other vertebrate
# 9 Invertebrate
# 12 Primate
# 14 Other mammalian
# 16 Rodent

# Connect to SQL database in restez path
restez_connect()

# Create the SQL database SLOW BE PATIENT
# set minimum sequence length to 500 bp
db_create(db_type = "nucleotide", min_length = 500, max_length = NULL,
          alt_restez_path = NULL)

 # Disconnect
restez_disconnect()













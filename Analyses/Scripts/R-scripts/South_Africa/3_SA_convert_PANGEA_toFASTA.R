# script to convert all sequence data into FASTA format
# this sequence will be later used to create a database to blast sequences
# locally

# I will exclude all sequences from South Africa and in which country is unknown

library(phylotools)

# Get the latest extract of the PANGEA data
pangea_data <- read.csv("~/Box Sync/pangea2018/data/PANGEA_Extract_Volz_2019-04-16.csv",
                        na.strings=c("","NA"))
# subset PANGEA data and get only sequences in which location is not South Africa
excSA <- subset(pangea_data, main_cohort_id != "South Africa" & geo_country != "SouthAfrica")

# get a datafrafe containing only sequences and pangea_id
new_data <- excSA[c("pangea_id", "shiver_consensus")]

# function requires to change name of columns
colnames(new_data) <- c("seq.name", "seq.text")

dat2fasta(new_data, outfile = "pangea_excSA.fasta")

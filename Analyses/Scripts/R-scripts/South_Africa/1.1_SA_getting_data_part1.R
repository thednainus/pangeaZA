# Get FASTA sequences for South Africa data using the latest PANGEA extract
# This script will get the full genomes for non-duplicated sequences

library(phylotools)
library(stringr)

################################################################################
### PART 1
################################################################################

# Read csv file with the lates PANGEA extract
pangea_data <- read.csv("~/Box Sync/pangea/data/PANGEA_Extract_Volz_2019-04-16.csv",
                        na.strings=c("","NA"))

# Subset data and get only sequences from South Africa
sa <- subset(pangea_data, main_cohort_id == "South Africa" | geo_country == "SouthAfrica")
# count number of base pairs and exclude gaps
sa["consensus_length"] <- nchar(sa$shiver_consensus) - str_count(sa$shiver_consensus, "-")


# Remove duplicated sequences using field participant_id
# I used the solution as described here
# https://stackoverflow.com/questions/25962909/remove-duplicates-based-on-2nd-column-condition
sa_dups <- sa[with(sa, ave(consensus_length, participant_id, FUN = max) == consensus_length),]
# However, sa_dups still contains some duplicates in case that consensus_length
# was identical between the duplicates
# to remove the rest of the duplicates, I do the below
sa_noDups <- sa_dups[!duplicated(sa_dups$participant_id),]


# save data as FASTA format
SA_all_data <- sa_noDups[c("pangea_id", "shiver_consensus")]
colnames(SA_all_data) <- c("seq.name", "seq.text")
dat2fasta(SA_all_data, outfile = "SA_all_data.fasta")

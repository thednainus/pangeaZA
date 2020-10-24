# reads results from Python script (blast_source_sequence_pangea.py.
# This Python script selectes the best hit for each South Africa sequence to be
# used as source sequence in pyhlodynamic analyses

library(phylotools)

# get pangea data to be able to get sequences for best matches
pangea_data <- read.csv("~/Box Sync/pangea2018/data/PANGEA_Extract_Volz_2019-04-16.csv",
                        na.strings=c("","NA"))


################################################################################
# GAG GENE

# reads results for pol gene
gag_bestHit <- read.csv("~/Box Sync/my_R_packages/pangea/Analyses/Blast_source/South_Africa/bestMatchIDs/SA_gag_id_3bestMatches_pangea.csv",
                        header = FALSE)
colnames(gag_bestHit) <- "bestHit"

# get the subset of the pangea_data dataframe that matches the best hits for
# gag gene (as in dataframe pol_bestHit)
pangea_subset_gag <- pangea_data[match(gag_bestHit$bestHit,
                                       pangea_data$pangea_id), ]
pangea_subset_gag <- pangea_subset_gag[c("pangea_id", "geo_country", "shiver_consensus")]


# save these sequences (from pangea_gag) as FASTA so I can add them in the
# alignment for the pol gene and use as source (src) sequences

# get a datafrafe containing only sequences and pangea_id

gag <- pangea_subset_gag[c(1,3)]

# function requires to change name of columns
colnames(gag) <- c("seq.name", "seq.text")
dat2fasta(gag, outfile = "SA_src_pangea_gag_3bestHits.fasta")


################################################################################
# POL GENE

# reads results for pol gene
pol_bestHit <- read.csv("~/Box Sync/my_R_packages/pangea/Analyses/Blast_source/South_Africa/bestMatchIDs/SA_pol_id_3bestMatches_pangea.csv",
                        header = FALSE)
colnames(pol_bestHit) <- "bestHit"

# get the subset of the pangea_data dataframe that matches the best hits for
# pol gene (as in dataframe pol_bestHit)
pangea_subset_pol <- pangea_data[match(pol_bestHit$bestHit,
                                       pangea_data$pangea_id), ]
pangea_subset_pol <- pangea_subset_pol[c("pangea_id", "geo_country", "shiver_consensus")]


# save these sequences (from pangea_pol) as FASTA so I can add them in the
# alignment for the pol gene and use as source (src) sequences

# get a datafrafe containing only sequences and pangea_id

pol <- pangea_subset_pol[c(1,3)]

# function requires to change name of columns
colnames(pol) <- c("seq.name", "seq.text")

dat2fasta(pol, outfile = "SA_src_pangea_pol_3bestHits.fasta")


################################################################################
# ENV GENE (when blasting the alignment rather than the sequences as it is.
# This is because too many best hits were observed for the sequences as it is.
# Maybe because of the hypervariable regions?
# In the alignment I removed these hypervariable regions.)

# reads results for env gene
env_ali_bestHit <- read.csv("~/Box Sync/my_R_packages/pangea/Analyses/Blast_source/South_Africa/bestMatchIDs/SA_env_ali_id_3bestMatches_pangea.csv",
                            header = FALSE)
colnames(env_ali_bestHit) <- "bestHit"

# get the subset of the pangea_data dataframe that matches the best hits for
# env gene (as in dataframe env_bestHit)
pangea_subset_env_ali <- pangea_data[match(env_ali_bestHit$bestHit,
                                           pangea_data$pangea_id), ]
pangea_subset_env_ali <- pangea_subset_env_ali[c("pangea_id", "geo_country", "shiver_consensus")]


# save these sequences (from pangea_pol) as FASTA so I can add them in the
# alignment for the pol gene and use as source (src) sequences

# get a datafrafe containing only sequences and pangea_id
env_ali <- pangea_subset_env_ali[c(1,3)]

# function requires to change name of columns
colnames(env_ali) <- c("seq.name", "seq.text")
dat2fasta(env_ali, outfile = "SA_src_pangea_env_ali_3bestHits.fasta")

save(gag, pol, env_ali, file = "SA_pangea_cgr.rda")

# Script to get data for env gene to be used with treedater
library(lubridate)
library(senegalHIVmodel)
library(ape)
library(pangea)

# Get data for South African sequences
source("Analyses/Scripts/R-scripts/South_Africa/get_SA_data.R")

#### convert dates to decimal for CGR (close global reference) GenBank sequences
# read csv files

env_cgr <- read.csv(system.file("data/SA/CGR_genbank_results/src_genbank_env.csv",
                                package = 'pangea'))

# using function from package senegalHIVmodel to convert dates to decimal
###########
# env gene
env_cgr_d <-data_format(env_cgr, "CGR")
env_cgr_df <- env_cgr_d[[1]]
env_cgr_ms <- env_cgr_d[[2]] # ms = missing sample

#Dataframe to be used in treedater for the missing sample date
env_ms_df <- env_cgr_ms[c("lower", "upper")]
row.names(env_ms_df) <- env_cgr_ms[,1]


# read rooted plylogenetic trees for env
# note that trees were rooted using FigTree version 1.4.3
#env_tr <- read.tree(system.file( "data/SA/recombinants_included/env_ML.tree", package = 'pangea' ))

# drop outgroup sequence
#env_tr2 <- drop.tip(env_tr, tip = c("K03454_env", "AY371157_env"))

# read rooted plylogenetic tree for env
# Here the tree is unrooted and I will root it using R before droping
# outgroup sequences
env_tr <- read.tree( system.file( "data/SA/no_recombinants/env_SA_CGR_ML.tre",
                                  package = 'pangea' ))

# Root phylogenetic tree and drop outgroup sequences
env_tr2 <- root_and_drop_tips(tree = env_tr,
                              outgroups = c("K03454_env", "AY371157_env"))


# get sequences that is in the phylogenetic tree for the SA sequences
# that best match is with GenBank sequences
env_dates <- env_cgr_df[env_cgr_df$Sequence_name %in% env_tr2$tip.label,]


# get sequences that is in the phylogenetic tree for pangea CGR
# the lapply function will simply remove the "_CGR" from some phylogenetic tip
# label to match sequences in the pangea dataset.
pangea_env <- pangea_data[pangea_data$pangea_id %in%
                            lapply(env_tr2$tip.label,
                                   function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1]),]
#remove sequences from SA because here I am interested in only CGR sequences
pangea_env_cgr <- subset(pangea_env, geo_country != "South Africa")

# get sample_date for each South African sequence
pangea_env_cgr_dates <- pangea_env_cgr[c(1,5,6,11)]

# check if there is any NA in the sample dates
sum(is.na(pangea_env_cgr_dates$sample_date)) # no NA in sample dates

# convert dates to year-month-day format
pangea_env_cgr_dates["cgr_ymd"] <- dmy(pangea_env_cgr_dates$sample_date)

# convert dates to decimal format (because treedater used dates in decimal
# format)
pangea_env_cgr_dates["decimal"] <- decimal_date(pangea_env_cgr_dates$cgr_ymd)
pangea_env_cgr_dates$decimal <- as.factor(pangea_env_cgr_dates$decimal)

# get all SA sequences in the phylogenetic tree
sa_env <- sa_dates[sa_dates$pangea_id %in% env_tr2$tip.label,]

# merge all sequences in the phylogenetic tree with the CGR genbank sequences
# note that some pangea data will also be used as CGR
env_sa_all <-sa_env[c(1,6)]
names(env_sa_all)[1] <- "Sequence_name"

env_pangea_cgr <- pangea_env_cgr_dates[c(1,6)]
names(env_pangea_cgr)[1] <- "Sequence_name"
env_pangea_cgr$Sequence_name <- paste(env_pangea_cgr$Sequence_name, "CGR", sep = "_")

env_genbank_cgr <- env_dates[c(1,19)]

# merge all data
all_env_dates <- rbind(env_sa_all, env_pangea_cgr, env_genbank_cgr)
all_env_dates$Sequence_name <- as.character(all_env_dates$Sequence_name)
all_env_dates$Sequence_name <- as.factor(all_env_dates$Sequence_name)

all_env_dates$decimal <- as.character(all_env_dates$decimal)
all_env_dates$decimal <- as.numeric(all_env_dates$decimal)



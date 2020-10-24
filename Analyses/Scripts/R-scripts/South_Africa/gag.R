# Script to get data for gag gene to be used with treedater
library(lubridate)
library(senegalHIVmodel)
library(ape)
library(pangea)

# Get data for South African sequences
source("Analyses/Scripts/get_SA_data.R")

#### convert dates to decimal for CGR (close global reference) GenBank sequences
# read csv files
gag_cgr <- read.csv(system.file("data/SA/CGR_genbank_results/src_genbank_gag.csv",
                                package = 'pangea'))


# using function from package senegalHIVmodel to convert dates to decimal
###########
# gag gene
gag_cgr_d <-data_format(gag_cgr, "CGR")
gag_cgr_df <- gag_cgr_d[[1]]
gag_cgr_ms <- gag_cgr_d[[2]] # ms = missing sample

#Dataframe to be used in treedater for the missing sample date
gag_ms_df <- gag_cgr_ms[c("lower", "upper")]
row.names(gag_ms_df) <- gag_cgr_ms[,1]


# read rooted plylogenetic tree for gag
# Here the tree is unrooted and I will root it using R before droping
# outgroup sequences
gag_tr <- read.tree( system.file( "data/SA/no_recombinants/gag_SA_CGR_ML.tre",
                                  package = 'pangea' ))

# Root phylogenetic tree and drop outgroup sequences
gag_tr2 <- root_and_drop_tips(tree = gag_tr,
                              outgroups = c("K03454_gag", "AY371157_gag"))


# get sequences that is in the phylogenetic tree for the SA sequences
# that best match is with GenBank sequences
gag_dates <- gag_cgr_df[gag_cgr_df$Sequence_name %in% gag_tr2$tip.label,]


# get sequences that is in the phylogenetic tree for pangea CGR
# the lapply function will simply remove the "_CGR" from some phylogenetic tip
# label to match sequences in the pangea dataset.
pangea_gag <- pangea_data[pangea_data$pangea_id %in%
                            lapply(gag_tr2$tip.label,
                                   function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1]),]
#remove sequences from SA because here I am interested in only CGR sequences
pangea_gag_cgr <- subset(pangea_gag, geo_country != "South Africa")

# get sample_date for each South African sequence
pangea_gag_cgr_dates <- pangea_gag_cgr[c(1,5,6,11)]

# check if there is any NA in the sample dates
sum(is.na(pangea_gag_cgr_dates$sample_date)) # no NA in sample dates

# convert dates to year-month-day format
pangea_gag_cgr_dates["cgr_ymd"] <- dmy(pangea_gag_cgr_dates$sample_date)

# convert dates to decimal format (because treedater uses dates in decimal
# format)
pangea_gag_cgr_dates["decimal"] <- decimal_date(pangea_gag_cgr_dates$cgr_ymd)
pangea_gag_cgr_dates$decimal <- as.factor(pangea_gag_cgr_dates$decimal)

# get all SA sequences in the phylogenetic tree
sa_gag <- sa_dates[sa_dates$pangea_id %in% gag_tr2$tip.label,]

# merge all sequences in the phylogenetic tree with the CGR genbank sequences
# note that some pangea data will also be used as CGR
gag_sa_all <-sa_gag[c(1,6)]
names(gag_sa_all)[1] <- "Sequence_name"
gag_pangea_cgr <- pangea_gag_cgr_dates[c(1,6)]
names(gag_pangea_cgr)[1] <- "Sequence_name"
gag_pangea_cgr$Sequence_name <- paste(gag_pangea_cgr$Sequence_name, "CGR", sep = "_")
gag_genbank_cgr <- gag_dates[c(1,19)]

# merge all data
all_gag_dates <- rbind(gag_sa_all, gag_pangea_cgr, gag_genbank_cgr)
all_gag_dates$Sequence_name <- as.character(all_gag_dates$Sequence_name)
all_gag_dates$Sequence_name <- as.factor(all_gag_dates$Sequence_name)

all_gag_dates$decimal <- as.character(all_gag_dates$decimal)
all_gag_dates$decimal <- as.numeric(all_gag_dates$decimal)

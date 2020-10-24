# Script to get data for pol gene to be used with treedater
library(lubridate)
library(senegalHIVmodel)
library(ape)
library(pangea)

# Get data for South African sequences
source("Analyses/Scripts/R-scripts/South_Africa/get_SA_data2.R")

#### convert dates to decimal for CGR (close global reference) GenBank sequences
# read csv files

pol_cgr <- read.csv(system.file("data/SA/CGR_genbank_results/SA_src_metadata_pol_new.csv",
                                package = 'pangea'))

# using function from package senegalHIVmodel to convert dates to decimal
###########
# pol gene
pol_cgr_d <-data_format(pol_cgr, "CGR")
pol_cgr_df <- pol_cgr_d[[1]]
pol_cgr_ms <- pol_cgr_d[[2]] # ms = missing sample

#Dataframe to be used in treedater for the missing sample date
pol_ms_df <- pol_cgr_ms[c("lower", "upper")]
row.names(pol_ms_df) <- pol_cgr_ms[,1]

# read plylogenetic tree for pol
# Here the tree is unrooted and I will root it using R before droping
# outgroup sequences
pol_tr <- read.tree( system.file( "data/SA/no_recombinants/pol_SA_CGR_ML_new.tre",
                                  package = 'pangea' ))

# Root phylogenetic tree and drop outgroup sequences
pol_tr2 <- root_and_drop_tips(tree = pol_tr,
                              outgroups = c("K03454_pol", "AY371157_pol"))

# dropping tips that should not be in the tree because of new database
# provided on 16/04/2019
tips_to_drop <- sa$pangea_id[!sa$pangea_id %in% sa2$pangea_id]
tips_to_drop <- as.character(tips_to_drop)

pol_tr3 <- ape::drop.tip(pol_tr2, tip = tips_to_drop)

# get sequences that is in the phylogenetic tree for the SA sequences
# that best match is with GenBank sequences
pol_dates <- pol_cgr_df[pol_cgr_df$Sequence_name %in% pol_tr3$tip.label,]


# get sequences that is in the phylogenetic tree for pangea CGR
# the lapply function will simply remove the "_CGR" from some phylogenetic tip
# label to match sequences in the pangea dataset.
pangea_pol <- pangea_data3[pangea_data3$pangea_id %in%
                            lapply(pol_tr3$tip.label,
                                   function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1]),]

#remove sequences from SA because here I am interested in only CGR sequences
pangea_pol_cgr <- subset(pangea_pol, geo_country != "South Africa")

# get sample_date for each cgr sequence
pangea_pol_cgr_dates <- pangea_pol_cgr[c("pangea_id","age","gender","sample_date")]

# check if there is any NA in the sample dates
sum(is.na(pangea_pol_cgr_dates$sample_date)) # 1 NA in sample dates
# get which sequence have sample_date = NA
pangea_cgr_na <- subset(pangea_pol_cgr, is.na(sample_date))[[1]]

# drop from tree the sequence from CGR that has sample date as NA
pol_tr3 <- ape::drop.tip(pol_tr3, tip = paste(pangea_cgr_na, "CGR", sep = "_"))

# remove pangea_id == pangea_cgr_na from pangea_pol_cgr_dates
pangea_pol_cgr_dates <- subset(pangea_pol_cgr_dates, pangea_id != pangea_cgr_na)


# convert dates to year-month-day format
pangea_pol_cgr_dates["cgr_ymd"] <- dmy(pangea_pol_cgr_dates$sample_date)

# convert dates to decimal format (because treedater used dates in decimal
# format)
pangea_pol_cgr_dates["decimal"] <- decimal_date(pangea_pol_cgr_dates$cgr_ymd)
pangea_pol_cgr_dates$decimal <- as.factor(pangea_pol_cgr_dates$decimal)

# get all SA sequences in the phylogenetic tree
sa_pol <- sa_dates[sa_dates$pangea_id %in% pol_tr3$tip.label,]

# merge all sequences in the phylogenetic tree with the CGR genbank sequences
# note that some pangea data will also be used as CGR
pol_sa_all <-sa_pol[c("pangea_id","decimal")]
names(pol_sa_all)[1] <- "Sequence_name"

pol_pangea_cgr <- pangea_pol_cgr_dates[c("pangea_id", "decimal")]
names(pol_pangea_cgr)[1] <- "Sequence_name"
pol_pangea_cgr$Sequence_name <- paste(pol_pangea_cgr$Sequence_name, "CGR", sep = "_")

pol_genbank_cgr <- pol_dates[c("Sequence_name","decimal")]

# merge all data
all_pol_dates <- rbind(pol_sa_all, pol_pangea_cgr, pol_genbank_cgr)
all_pol_dates$Sequence_name <- as.character(all_pol_dates$Sequence_name)
all_pol_dates$Sequence_name <- as.factor(all_pol_dates$Sequence_name)

all_pol_dates$decimal <- as.character(all_pol_dates$decimal)
all_pol_dates$decimal <- as.numeric(all_pol_dates$decimal)

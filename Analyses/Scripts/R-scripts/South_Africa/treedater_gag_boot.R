# Function to read csv file containing metadata (for HIV sampling dates)
# and format it in a way to be used as input data for treedater

library(treedater)
library(lubridate)
library(phytools)


source("Analyses/Scripts/R-scripts/gag.R")

# Rooting bootstrap trees using the outgroups and removing outgroups from
# phylogenetic trees
outgroups <- c("K03454_gag", "AY371157_gag")

gag_boot <- read.tree(system.file("data/SA/no_recombinants/gag_SA_CGR_100bootstrap.tre",
                                  package = "pangea"))
gag_boot_new <- lapply(unclass(gag_boot), root_and_drop_tips, outgroups = outgroups)

# drop outlier tips from bootstrapped trees
# These outliers were identified when using treedater on the best ML tree
ot0 <- readRDS(system.file("data/SA/no_recombinants/gag_outliers.RDS",
                           package = "pangea"))
gag_boot_new2 <- lapply(gag_boot_new, drop.tip, as.character( ot0$taxon[ ot0$q < .05] ))

class(gag_boot_new2)<-"multiPhylo"

################# TREEDATER ON BOOTSTRAPPED TREES #############################
# load treedater tree for gag gene
load(system.file("data/SA/no_recombinants/gag_treedater.rda",
                 package = "pangea"))



boot_gag_dtr <- boot(td = dtr_gag2, tres = gag_boot_new2, ncpu = 2, quiet = FALSE,
                     overrideTempConstraint = FALSE)

saveRDS(boot_gag_dtr, file = "gag_bootstrap_treedater.RDS")

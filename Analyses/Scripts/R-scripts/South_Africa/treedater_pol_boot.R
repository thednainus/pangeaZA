# Function to read csv file containing metadata (for HIV sampling dates)
# and format it in a way to be used as input data for treedater

library(treedater)
library(lubridate)
library(phytools)


source("Analyses/Scripts/R-scripts/South_Africa/pol.R")

# Rooting bootstrap trees using the outgroups and removing outgroups from
# phylogenetic trees
outgroups <- c("K03454_pol", "AY371157_pol")

pol_boot <- read.tree(system.file("data/SA/no_recombinants/pol_SA_CGR_100bootstrap.tre",
                                  package = "pangea"))
pol_boot_new <- lapply(unclass(pol_boot), root_and_drop_tips, outgroups = outgroups)

# drop outlier tips from bootstrapped trees
# These outliers were identified when using treedater on the best ML tree
ot0 <- readRDS(system.file("data/SA/no_recombinants/pol_outliers.RDS",
                           package = "pangea"))
pol_boot_new2 <- lapply(pol_boot_new, drop.tip, as.character( ot0$taxon[ ot0$q < .05] ))

class(pol_boot_new2)<-"multiPhylo"

################# TREEDATER ON BOOTSTRAPPED TREES #############################
# load treedater tree for pol gene
load(system.file("data/SA/no_recombinants/pol_treedater.rda",
                           package = "pangea"))



boot_pol_dtr <- boot(td = dtr_pol2, tres = pol_boot_new2, ncpu = 2, quiet = FALSE,
                     overrideTempConstraint = FALSE)

saveRDS(boot_pol_dtr, file = "pol_bootstrap_treedater.RDS")

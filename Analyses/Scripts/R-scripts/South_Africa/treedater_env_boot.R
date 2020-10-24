# Function to read csv file containing metadata (for HIV sampling dates)
# and format it in a way to be used as input data for treedater

library(treedater)
library(lubridate)
library(phytools)


source("Analyses/Scripts/env.R")

# Rooting bootstrap trees using the outgroups and removing outgroups from
# phylogenetic trees
outgroups <- c("K03454_env", "AY371157_env")

env_boot <- read.tree(system.file("data/SA/no_recombinants/env_SA_CGR_100bootstrap.tre",
                                  package = "pangea"))
env_boot_new <- lapply(unclass(env_boot), root_and_drop_tips, outgroups = outgroups)

# drop outlier tips from bootstrapped trees
# These outliers were identified when using treedater on the best ML tree
ot0 <- readRDS(system.file("data/SA/no_recombinants/env_outliers.RDS",
                           package = "pangea"))
env_boot_new2 <- lapply(env_boot_new, drop.tip, as.character( ot0$taxon[ ot0$q < .05] ))

class(env_boot_new2)<-"multiPhylo"

################# TREEDATER ON BOOTSTRAPPED TREES #############################
# load treedater tree for env gene
load(system.file("data/SA/no_recombinants/env_treedater.rda",
                 package = "pangea"))



boot_env_dtr <- boot(td = dtr_env2, tres = env_boot_new2, ncpu = 2, quiet = FALSE,
                     overrideTempConstraint = FALSE)

saveRDS(boot_env_dtr, file = "env_bootstrap_treedater.RDS")

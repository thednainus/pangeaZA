# script to convert branch length of phylogenetic trees into units of
# calendar time using treedater (https://github.com/emvolz/treedater)
library(treedater)
library(lubridate)
library(ape)
library(pangea)
library(senegalHIVmodel)

#PACKAGEDIR <- './'

source("Analyses/Scripts/env.R")

########################### TREEDATER ##########################################

env.seqlen <- 2403 # the length of the HIV sequences

env.sts <- all_env_dates$decimal
names(env.sts) <- all_env_dates$Sequence_name

(dtr_env <- dater(env_tr2, env.sts, env.seqlen, estimateSampleTimes = env_ms_df, ncpu = 2 ))

# check outliers
outlierTips( dtr_env )  -> ot0
env_tr3 <- drop.tip( env_tr2, as.character( ot0$taxon[ ot0$q < .05] ) )
dtr_env2 <- dater(env_tr3, env.sts, env.seqlen, estimateSampleTimes = env_ms_df, omega0 = dtr_env$meanRate , ncpu = 20 )

# NOTE AY894994_CGR pushes tmrca back about 30 years
invisible('
> dtr_env2

Phylogenetic tree with 1399 tips and 1398 internal nodes.

Tip labels:
	AY894994_CGR, KX907421_CGR, PG14-ZA100188, PG14-ZA100939, KX228818_CGR, PG14-ZA102850, ...

Rooted; includes branch lengths.

 Time of common ancestor
1930.58343263827

 Time to common ancestor (before most recent sample)
84.8987591425523

 Mean substitution rate
0.00501063607732793

 Strict or relaxed clock
relaxed

 Coefficient of variation of rates
0.170300308456751

> max( node.depth.edgelength( drop.tip( dtr_env2, "AY894994_CGR" ) ))
[1] 52.87736
> max( na.omit(env.sts) ) - max( node.depth.edgelength( drop.tip( dtr_env2, "AY894994_CGR" ) ))
[1] 1962.605

> rootToTipRegressionPlot( dtr_env2  )
Root-to-tip mean rate: 0.00621642975459415
Root-to-tip p value: 4.24249245872444e-56
Root-to-tip R squared (variance explained): 0.163323516621434

# quick nonparametric phylodynamic analysis
library(skygrowth)
tr <- dtr_env2
class(tr) <- "phylo"
tr <- drop.tip( tr, tr$tip.label[ grepl("CGR", tr$tip.label ) ] )
f <- skygrowth.map( tr )
f$time <- f$time + max( na.omit( env.sts  ) )
plot( f, logy=FALSE ) + ggplot2::xlim( c(1975, 2015 ) )

')

save(dtr_env2, file = "env_treedater.rda")

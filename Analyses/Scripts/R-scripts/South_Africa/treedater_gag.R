# script to convert branch length of phylogenetic trees into units of
# calendar time using treedater (https://github.com/emvolz/treedater)
library(treedater)
library(lubridate)
library(ape)
library(pangea)
library(senegalHIVmodel)
library(phytools)

#PACKAGEDIR <- '/home/erik/git/pangea/'

source("Analyses/Scripts/R-scrits/gag.R")

########################### TREEDATER ##########################################
## gag gene

gag.seqlen <- 1476 # the length of the HIV sequences

gag.sts <- all_gag_dates$decimal
names(gag.sts) <- all_gag_dates$Sequence_name

(dtr_gag <- dater(gag_tr2, gag.sts, gag.seqlen, estimateSampleTimes = gag_ms_df, ncpu = 1 ))

ot0 <- outlierTips( dtr_gag )
saveRDS(ot0, "gag_outliers.RDS")
invisible('                      taxon            q            p    loglik       rates
PG14-ZA100135 PG14-ZA100135 0.0005600194 3.827884e-07 -9.857344 0.003660293
              branch.length
PG14-ZA100135             0
')

gag_tr3 <- drop.tip( gag_tr2, as.character( ot0$taxon[ ot0$q < .05] ) )
dtr_gag2 <- dater(gag_tr3, gag.sts, gag.seqlen, estimateSampleTimes = gag_ms_df, omega0 = dtr_gag$meanRate , ncpu = 2 )


save(dtr_gag2, file = "gag_treedater.rda")

# Results when recombinant sequences were presented on the the phylogenetic tree
invisible('
> dtr_gag

Phylogenetic tree with 1463 tips and 1462 internal nodes.

Tip labels:
	PG14-ZA101878, PG14-ZA102626, PG14-ZA101151, PG14-ZA100939, PG14-ZA100031, PG14-ZA100412, ...

Rooted; includes branch lengths.

 Time of common ancestor
1966.68136929536

 Time to common ancestor (before most recent sample)
48.8008224854648

 Mean substitution rate
0.00347443839936578

 Strict or relaxed clock
relaxed

 Coefficient of variation of rates
0.206996952636277

> dtr_gag2

Phylogenetic tree with 1462 tips and 1461 internal nodes.

Tip labels:
	PG14-ZA101878, PG14-ZA102626, PG14-ZA101151, PG14-ZA100939, PG14-ZA100031, PG14-ZA100412, ...

Rooted; includes branch lengths.

 Time of common ancestor
1964.49816401083

 Time to common ancestor (before most recent sample)
50.9840277699921

 Mean substitution rate
0.00331264715093059

 Strict or relaxed clock
relaxed

 Coefficient of variation of rates
0.206008574017109


Root-to-tip mean rate: 0.00375530321595699
Root-to-tip p value: 1.95727971652205e-47
Root-to-tip R squared (variance explained): 0.133593718742268



# quick nonparametric phylodynamic analysis
library(skygrowth)
tr <- dtr_gag2
class(tr) <- "phylo"
tr <- drop.tip( tr, tr$tip.label[ grepl("CGR", tr$tip.label ) ] )
f <- skygrowth.map( tr )
f$time <- f$time + max( na.omit( gag.sts  ) )
plot( f, logy=FALSE ) + ggplot2::xlim( c(1975, 2015 ) )

')


# Results after re-estimating tree without recombinant sequences
invisible('
> dtr_gag

          Phylogenetic tree with 1461 tips and 1460 internal nodes.

          Tip labels:
          PG14-ZA101299, PG14-ZA101994, PG14-ZA102074, PG14-ZA102736, PG14-ZA102508, PG14-ZA102576, ...

          Rooted; includes branch lengths.

          Time of common ancestor
          1975.07202443045

          Time to common ancestor (before most recent sample)
          40.410167350366

          Mean substitution rate
          0.00396027476375376

          Strict or relaxed clock
          relaxed

          Coefficient of variation of rates
          0.205014258671369

          > dtr_gag2

          Phylogenetic tree with 1459 tips and 1458 internal nodes.

          Tip labels:
          PG14-ZA101299, PG14-ZA101994, PG14-ZA102074, PG14-ZA102736, PG14-ZA102508, PG14-ZA102576, ...

          Rooted; includes branch lengths.

          Time of common ancestor
          1971.87680153703

          Time to common ancestor (before most recent sample)
          43.6053902437902

          Mean substitution rate
          0.00360203793218317

          Strict or relaxed clock
          relaxed

          Coefficient of variation of rates
          0.20232774619547
          ')





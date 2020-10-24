# script to convert branch length of phylogenetic trees into units of
# calendar time using treedater (https://github.com/emvolz/treedater)
library(treedater)
library(lubridate)
library(ape)
library(pangea)
library(senegalHIVmodel)

#PACKAGEDIR <- './'
#PACKAGEDIR <- '/home/erik/git/pangea/'


source("Analyses/Scripts/R-scripts/South_Africa/pol.R")

########################### TREEDATER ##########################################
## pol gene

pol.seqlen <- 2733 # the length of the HIV sequences

pol.sts <- all_pol_dates$decimal
names(pol.sts) <- all_pol_dates$Sequence_name

(dtr_pol <- dater(pol_tr3, pol.sts, pol.seqlen, estimateSampleTimes = pol_ms_df))


ot0 <- outlierTips( dtr_pol )
saveRDS(ot0, file = "pol_outliers_new.RDS")
# Results before recombinats were removed from alignments
invisible('
                              taxon            q            p     loglik       rates
          PG17-ZA000129         PG17-ZA000129 0.0000195881 2.183734e-08 -13.227179 0.002895664
          PG17-ZA001862         PG17-ZA001862 0.0003806319 8.486777e-07 -14.672734 0.001368103
          PG17-ZA000187         PG17-ZA000187 0.0008006149 2.677642e-06 -13.453017 0.004523279
          PG17-ZA001882         PG17-ZA001882 0.0008126372 3.623800e-06  -8.365138 0.002869798
          PG16-BW003356_CGR PG16-BW003356_CGR 0.0010862973 7.266203e-06 -12.820431 0.004691773
          PG17-ZA000124         PG17-ZA000124 0.0010862973 6.950811e-06 -12.792375 0.001278753
          PG14-ZA101031         PG14-ZA101031 0.0095769797 8.541342e-05 -10.429149 0.001396785
          PG14-ZA101926         PG14-ZA101926 0.0095769797 7.735449e-05 -11.517531 0.004925634
          PG15-BW001023_CGR PG15-BW001023_CGR 0.0212510316 2.132211e-04 -11.059048 0.001154447
          PG16-BW003355_CGR PG16-BW003355_CGR 0.0222675442 2.482446e-04  -8.994243 0.001749637
          branch.length
          PG17-ZA000129          0.000000
          PG17-ZA001862          2.778695
          PG17-ZA000187          1.311041
          PG17-ZA001882          0.000000
          PG16-BW003356_CGR      1.960464
          PG17-ZA000124          3.931238
          PG14-ZA101031          3.942488
          PG14-ZA101926          6.013721
          PG15-BW001023_CGR     28.817249
          PG16-BW003355_CGR      1.490519
')


# Results after recombinats were removed from alignments
invisible('
                                    taxon            q            p      loglik     rates          branch.length
          PG16-UG000765_CGR PG16-UG000765_CGR 1.011506e-15 6.201751e-19 -43.09543 0.0007832606     46.403121
          PG14-ZA101511         PG14-ZA101511 8.199184e-10 1.005418e-12 -30.58772 0.0079267604     22.549977
          PG14-ZA101031         PG14-ZA101031 9.331916e-04 1.716478e-06 -14.17828 0.0011720793      4.601304
          PG14-ZA101946         PG14-ZA101946 1.195489e-02 2.931918e-05 -11.13042 0.0015378311      2.026108
          PG14-ZA102892         PG14-ZA102892 1.519223e-02 4.657337e-05 -13.00662 0.0050830809     18.496464
          ')

pol_tr3 <- drop.tip( pol_tr2, as.character( ot0$taxon[ ot0$q < .05] ) )
dtr_pol2 <- dater(pol_tr3, pol.sts, pol.seqlen, estimateSampleTimes = pol_ms_df, omega0 = dtr_pol$meanRate , ncpu = 20 )

# Results before removing recombinants from the phylogenetic tree
invisible("
> dtr_pol

NOTE: The p values for lineage clock rates show at least one outlying value after adjusting for multiple-testing.  This indicates a poor fit to the data for at least a portion of the phylogenetic tree. To visualize the distribution p-values, use *goodnessOfFitPlot*.

Phylogenetic tree with 1634 tips and 1633 internal nodes.

Tip labels:
	PG16-UG000765_CGR, PG14-ZA102949, PG14-ZA102793, PG14-ZA102545, PG14-ZA101803, AF110980_CGR, ...

Rooted; includes branch lengths.

 Time of common ancestor
1974.21083592473

 Time to common ancestor (before most recent sample)
41.7809448971875

 Mean substitution rate
0.00269362987713857

 Strict or relaxed clock
relaxed

 Coefficient of variation of rates
0.21364214068346

> dtr_pol2

Phylogenetic tree with 1629 tips and 1628 internal nodes.

Tip labels:
	PG14-ZA102949, PG14-ZA102793, PG14-ZA102545, PG14-ZA101803, AF110980_CGR, PG14-ZA102744, ...

Rooted; includes branch lengths.

 Time of common ancestor
1955.41106139181

 Time to common ancestor (before most recent sample)
60.0711303890102

 Mean substitution rate
0.0025421503592312

 Strict or relaxed clock
relaxed

 Coefficient of variation of rates
0.185876129149531


> rootToTipRegressionPlot( dtr_pol2 )
Root-to-tip mean rate: 0.00260497363018862
Root-to-tip p value: 9.22430828542404e-33
Root-to-tip R squared (variance explained): 0.0836605015758213

# quick nonparametric phylodynamic analysis
library(skygrowth)
tr <- dtr_pol2
class(tr) <- 'phylo'
tr <- drop.tip( tr, tr$tip.label[ grepl('CGR', tr$tip.label ) ] )
f <- skygrowth.map( tr )
f$time <- f$time + max( na.omit( pol.sts  ) )
plot( f, logy=FALSE ) + ggplot2::xlim( c(1975, 2015 ) )

")

save(dtr_pol2, file = "pol_treedater.rda")

# Results after removing recombinants from phylogenetic tree and re-estimating
# tree
invisible("
> dtr_pol

          NOTE: The p values for lineage clock rates show at least one outlying value after adjusting for multiple-testing.  This indicates a poor fit to the data for at least a portion of the phylogenetic tree. To visualize the distribution p-values, use *goodnessOfFitPlot*.

          Phylogenetic tree with 1631 tips and 1630 internal nodes.

          Tip labels:
          PG14-ZA101301, PG14-ZA101938, PG14-ZA102536, PG14-ZA101187, PG14-ZA101685, PG14-ZA101103, ...

          Rooted; includes branch lengths.

          Time of common ancestor
          1969.58866024213

          Time to common ancestor (before most recent sample)
          46.4031205797862

          Mean substitution rate
          0.00250792723311467

          Strict or relaxed clock
          relaxed

          Coefficient of variation of rates
          0.205945296336946


          > dtr_pol2

          Phylogenetic tree with 1626 tips and 1625 internal nodes.

          Tip labels:
          PG14-ZA101301, PG14-ZA101938, PG14-ZA102536, PG14-ZA101187, PG14-ZA101685, PG14-ZA101103, ...

          Rooted; includes branch lengths.

          Time of common ancestor
          1959.84275032316

          Time to common ancestor (before most recent sample)
          55.639441457655

          Mean substitution rate
          0.00230318341301805

          Strict or relaxed clock
          relaxed

          Coefficient of variation of rates
          0.191410994072377


          > rootToTipRegressionPlot( dtr_pol2 )
          Root-to-tip mean rate: 0.002576460303722
          Root-to-tip p value: 1.622507845158e-31
          Root-to-tip R squared (variance explained): 0.0805893434043493

")

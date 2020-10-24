#!/bin/bash
# script to perform some post modifications to the results of bootstrap analysis for ML trees
# ML = maximum likelihood trees
# Note that you may want to change some of the file names

# CHANGE to the directory of interest
cd //Your/Path/Here/Analyses/Trees/SA_and_src/3bestMatches/env/Bootstrap/results

#join raxml.bootstraps files and remove the rest
cat BOOT*.raxml.bootstraps > final.BOOT.SA.CGR.HKY_Gp12+3.raxml.bootstraps
rm BOOT*.bootstraps

#join raxml.log files and remove the rest
cat BOOT*.raxml.log > final.BOOT.SA.CGR.HKY_Gp12+3.raxml.log
rm BOOT*.raxml.log

#join startTrees files into a single file and remove the rest
#cat BOOT*.raxml.startTree > final.BOOT.SA.CGR.HKY_Gp12+3.raxml.startTree
#rm BOOT*.raxml.startTree



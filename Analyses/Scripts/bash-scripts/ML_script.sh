#!/bin/bash
# script to perform some post modifications to the results of analysis of best ML tree search
# get the best ML tree with the highest ML score
# ML = maximum likelihood tree
# Note that you may need to change file names in some lines

# CHANGE to the directory of interest
cd /Your/Path/Here/Analyses/Trees/SA_and_src/3bestMatches/env/ML/results

#deals with choosing the best likelihood tree
grep "Final LogLikelihood" *.log > loglikelihoods.txt
sort -t : -k 3 -g -r loglikelihoods.txt > sorted_loglikelihoods.txt

first_line=`head -n1 sorted_loglikelihoods.txt`

#gets first occurrence of $first_line using as separator the dot '.'
highest_ln=(${first_line//./ })

#gets best ML tree (with highest log likelihood value)
best_ML_tree=$highest_ln.raxml.bestTree
#mv $best_ML_tree final.$best_ML_tree

#remove the other bestTree files
#rm ML*.raxml.bestTree

#join bestTree files into a single file (mlTrees) and remove the rest
cat *.raxml.bestTree > final.ML.SA.CGR.HKY_Gp12+3.raxml.mlTrees
mv $best_ML_tree final.$best_ML_tree
rm ML*.raxml.bestTree

#save the bestModel for best tree in a single file
#and remove the rest
mv $highest_ln.raxml.bestModel final.ML.SA.CGR.HKY_Gp12+3.raxml.bestModel
rm ML*.raxml.bestModel

#join startTrees files into a single file and remove the rest
cat ML*.raxml.startTree > final.ML.SA.CGR.HKY_Gp12+3.raxml.startTree
rm ML*.raxml.startTree

#join log files into a single file and remove the rest
cat ML*.raxml.log > final.ML.SA.CGR.HKY_Gp12+3.raxml.log
rm ML*.raxml.log

#reduced files will be create in each run, but they are identical
#and I should keep just one of these files
mv ML_0_GTR+I+G.raxml.reduced.phy final.ML.SA.CGR.HKY_Gp12+3.raxml.reduced.phy
rm ML*.raxml.reduced.phy



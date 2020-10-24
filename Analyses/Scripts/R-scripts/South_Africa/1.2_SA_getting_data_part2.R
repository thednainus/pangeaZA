library(phylotools)

################################################################################
### PART 2
################################################################################

# Read files for each gene (pol, gag and env) and get only sequences
# legth >= 800 base pairs and save as fasta format
# Files SA_pol.csv, SA_gag.csv and SA_env.csv were obtained using
# blastn command line to align two sequences. Query sequence using a reference
# sequence for subtype C was used to separate each SA sequence into gag, pol
# and env genes. Then I got the sequence size and
# aligned bit from the subject sequence to create these csv files.
SA_pol <- read.csv("Analyses/Sequences/South_Africa/SA_pol.csv",
                   header = FALSE)

SA_gag <- read.csv("Analyses/Sequences/South_Africa/SA_gag.csv",
                   header = FALSE)

SA_env <- read.csv("Analyses/Sequences/South_Africa/SA_env.csv",
                   header = FALSE)

# subset data to get sequence lenght >= 800
SA_pol2 <- subset(SA_pol, V2 >= 800)
SA_gag2 <- subset(SA_gag, V2 >= 800)
SA_env2 <- subset(SA_env, V2 >= 800)


#save as fasta format
SA_pol3 <- SA_pol2[c(1,3)]
colnames(SA_pol3) <- c("seq.name", "seq.text")
dat2fasta(SA_pol3, outfile = "SA_pol.fasta")

SA_gag3 <- SA_gag2[c(1,3)]
colnames(SA_gag3) <- c("seq.name", "seq.text")
dat2fasta(SA_gag3, outfile = "SA_gag.fasta")

SA_env3 <- SA_env2[c(1,3)]
colnames(SA_env3) <- c("seq.name", "seq.text")
dat2fasta(SA_env3, outfile = "SA_env.fasta")


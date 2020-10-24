library(ape)

################################################################################
# Analysing results from using COMET HIV-1 (https://comet.lih.lu/) for
# subtyping sequences.
################################################################################

SA_gag_st <- read.delim("Analyses/Subtyping/SouthAfrica/SA_gag_COMET.csv")
SA_pol_st <- read.delim("Analyses/Subtyping/SouthAfrica/SA_pol_COMET.csv")
SA_env_st <- read.delim("Analyses/Subtyping/SouthAfrica/SA_env_COMET.csv")

# subset dataframes and check only sequences that are not from subtype C
# these sequences will be later removed from alignments
sSA_gag_st <- subset(SA_gag_st, subtype != "C")
sSA_pol_st <- subset(SA_pol_st, subtype != "C")
sSA_env_st <- subset(SA_env_st, subtype != "C")

################################################################################
# Remove sequences that are not from subtype C from the alignment

# Read gag, pol and env sequences
SA_gag_seq <- read.FASTA("Analyses/Sequences/South_Africa/byGene/SA_gag.fasta")
SA_pol_seq <- read.FASTA("Analyses/Sequences/South_Africa/byGene/SA_pol.fasta")
SA_env_seq <- read.FASTA("Analyses/Sequences/South_Africa/byGene/SA_env.fasta")


# get only sequences classified as subtype C
SA_gag_stC <- SA_gag_seq[!names(SA_gag_seq) %in% sSA_gag_st$name]
write.FASTA(SA_gag_stC, file = "SA_gag_C.fasta")

SA_pol_stC <- SA_pol_seq[!names(SA_pol_seq) %in% sSA_pol_st$name]
write.FASTA(SA_pol_stC, file = "SA_pol_C.fasta")

SA_env_stC <- SA_env_seq[!names(SA_env_seq) %in% sSA_env_st$name]
write.FASTA(SA_env_stC, file = "SA_env_C.fasta")

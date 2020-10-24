# Checking for drug resitance sites and mask them in a DNA alignment
# HIV sequences
library(ape)
library(big.phylo)
# mask drug resistance sites using the big.phylo package

seq.rm.drugresistance <- function (seq, outfile = NA)
{
  require(data.table)
  stopifnot(any(rownames(seq) == "HXB2"))
  load(system.file(package = "big.phylo", "AC_drugresistance_201508.rda"))
  load(system.file(package = "big.phylo", "refseq_hiv1_hxb2.rda"))
  hxb2 <- paste(hxb2[, HXB2.K03455], collapse = "")
  seq.hxb2 <- paste(as.character(seq["HXB2", ]), collapse = "")
  seq.hxb2.ng <- gsub("-+", "", seq.hxb2)
  stopifnot(nchar(seq.hxb2.ng) <= nchar(hxb2))
  seq.hxb2.st <- as.integer(regexpr(substr(seq.hxb2.ng, 1,
                                           20), hxb2))
  stopifnot(seq.hxb2.st > 0)
  stopifnot(nchar(seq.hxb2.ng) + seq.hxb2.st - 1L <= nchar(hxb2))
  stopifnot(seq.hxb2.ng == substr(hxb2, seq.hxb2.st, nchar(seq.hxb2.ng) +
                                    seq.hxb2.st - 1L))
  seq.hxb2.pos <- c()
  tmp <- gregexpr("-+", seq.hxb2)[[1]]
  if (tmp[1] > -1) {
    seq.hxb2.pos <- data.table(GP_ST = as.integer(tmp), GP_LEN = attr(tmp,
                                                                      "match.length"))
    seq.hxb2.pos <- seq.hxb2.pos[, list(GP_POS = seq.int(GP_ST,
                                                         length = GP_LEN)), by = "GP_ST"][, GP_POS]
  }
  seq.hxb2.pos <- data.table(HXB2INSEQ_POS = setdiff(seq_len(nchar(seq.hxb2)),
                                                     seq.hxb2.pos))
  seq.hxb2.pos[, `:=`(HXB2_POS, seq_len(nrow(seq.hxb2.pos)))]
  set(seq.hxb2.pos, NULL, "HXB2_POS", seq.hxb2.pos[, HXB2_POS +
                                                     seq.hxb2.st - 1L])
  stopifnot(seq.hxb2.pos[, tail(HXB2_POS, 1)] <= nchar(hxb2))
  setnames(dr, "HXB2.pos", "HXB2_POS")
  dr <- merge(dr, seq.hxb2.pos, by = "HXB2_POS")
  setnames(dr, "HXB2INSEQ_POS", "Alignment.nuc.pos")
  tmp <- seq.rm.drugresistance.internal(as.character(seq),
                                        dr, verbose = 1, rtn.DNAbin = 1)
  dr.info <- tmp$nodr.info
  nodr.seq <- tmp$nodr.seq
  if (!is.na(outfile))
    save(seq, nodr.seq, dr.info, file = outfile)
  list(nodr.info = dr.info, nodr.seq = nodr.seq)
}

seq.rm.drugresistance.internal<- function(char.matrix, dr, verbose=1, rtn.DNAbin=1)
{
  if(verbose)	cat(paste("\nchecking for potential drug resistance mutations, n=",nrow(dr)))
  nodr.info	<- dr[, {
    query.yes	<- seq.find(char.matrix, Alignment.nuc.pos, unlist(strsplit(unlist(Mutant.NTs),'')))
    #print (char.matrix)
    if(length(query.yes))
    {
      ans		<- rownames(char.matrix)[query.yes]
      #print( char.matrix[query.yes, seq.int(dr[i,Alignment.nuc.pos]-3, length.out=9)] ); stop()
    }
    if(!length(query.yes))
      ans		<- NA_character_
    list(TAXA=ans)
  }, by=c('HXB2_POS','DR.name','Gene.codon.number','Wild.type','Mutant.NTs','Alignment.nuc.pos')]
  nodr.info	<- subset(nodr.info, !is.na(TAXA))
  for(i in nodr.info[, sort(unique(Alignment.nuc.pos))])
  {
    tmp		<- subset(nodr.info, Alignment.nuc.pos==i)[, TAXA]
    stopifnot(length(tmp)>0)
    # I changed the nnn for ---, as my alignments should have gaps and not n
    # this change was made by me on 11/09/2017
    cat(paste('\nsetting at pos',i,'to --- for taxa, n=', length(tmp)))
    char.matrix[tmp,	seq.int(i, length.out=3) ]<- matrix("-", nrow=length(tmp), ncol=3)
  }
  if(rtn.DNAbin)
    char.matrix	<- as.DNAbin(char.matrix)
  list(nodr.seq=char.matrix, nodr.info=nodr.info)
}


# function to read a sequence alignment that must contain reference HIV seq HXB2
# data_dir = directory where alignment is saved and masked alignment will be saved

mask_seq <- function(ali_file_path, wd){
  # wording directory in which files should be saved
  setwd(wd)

  # read DNA alignment in R
  seq <- read.dna(ali_file_path, format="sequential")
  # get the file name from a path
  filename <- basename(ali_file_path)

  # run funtion in package big.phylo to get info on masked alignment and
  # drug resistance sites
  info.seq.dr <- seq.rm.drugresistance(seq)

  # write information on drug resistance sites that were masked
  out_filename <- gsub(".fasta", "", filename)
  write.csv(info.seq.dr$nodr.info,file=paste(out_filename, "_DRS.csv", sep=""))
  write.dna(info.seq.dr$nodr.seq, file=paste("masked", filename, sep="_"))

}

# wd = working directory
wd <- "~/Box Sync/my_R_packages/pangea/Analyses/Alignments/South_Africa/byCodon/drug_resistance"
input_files <- "~/Box Sync/my_R_packages/pangea/Analyses/Alignments/South_Africa/byCodon/drug_resistance/phylip_alignment"

#input_files <- "/Users/user/Box Sync/linkages/data/senegal/modified_by_Fabricia/DRS_masked/reference/"
# List all files that contain ".phy" in "input_files"
ali_files <- as.matrix(Sys.glob(file.path(input_files, "*.phy")))

sapply(ali_files[2], mask_seq, wd)

#mask_seq(ali_files[1], wd)

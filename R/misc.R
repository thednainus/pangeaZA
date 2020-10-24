#' Setup Model Equations
#'
#' Function to setup the components of the mathematical model. See \code{\link[phydynR]{build.demographic.process}}
#' for more details
#'
#' @param demes a character vector naming the demes of the mathematical model.
#' @param nondemes a character vector naming the non demes of the mathematical model.
#' @param rcpp if TRUE, the expressions are interpreted as C code using the Rcpp package.
#'
#' @return this function returns a list containing the empty components
#'  (represented by zeros) to build the mathematical model.
#'  These components are the birth, death, migrations, total number of demes and
#'  non-demes of the model.
#'  \itemize{
#'      \item Birth is a vector or matrix describing the model birth rates;
#'      \item Death is a vector describing the model death rates;
#'      \item Migration is a vector or matrix describing the model migration
#'                rates.
#'  }
#'
#'
#' @seealso \code{\link[phydynR]{build.demographic.process}}
#'
#' @export
#'
#' @examples
#' demes <- c("gpm", "gpf", "msm", "src")
#' eqns <- setup.model.equations(demes)
setup.model.equations <- function(demes, nondemes = NULL, rcpp = FALSE)
{
  # size of demes object
  m <- length(demes)
  # size of nondemes object
  mm <- length(nondemes)
  # birth matrix filled with zeros
  b <- matrix('0.', nrow = m, ncol = m)
  # migration matrix filled with zeros
  migs <- matrix('0.', nrow = m, ncol = m)

  rownames(b) = rownames(migs) = colnames(b) = colnames(migs) <- demes

  dths <- setNames(rep('0.', m ), demes)

  ndd <- setNames(rep('0.', mm ), nondemes)

  # return the components of the mathematical model as a list
  list(births = b, migs = migs, deaths = dths,
       nonDemeDynamics = ndd, m = m, mm = mm ,
       demes = demes, nondemes = nondemes)
}

#' Title Get fasta sequences from PANGEA dataframe
#'
#' Function to select sequences from specific country and save these DNA sequences
#' as FASTA.
#'
#' @param all_data Dataframe containing all data
#' @param country_name String for country name of interest
#' @param outfile_name String with the output filenama to save sequences
#' @param add_unknown value "yes" or "no". If "yes" include all genders (females,
#'  males and unknown), if "no" do not add sample in which gender is unknown.
#'
#' @return Fasta file for selected DNA sequences.
#' @export
#'
#' @examples #To Do.
get_sequences <- function(all_data, country_name, outfile_name, add_unknown){
  country_df <- subset(all_data, geo_country == country_name)

  if(add_unknown == "no"){
    country_df <- subset(country_df, gender != "Unknown")
  }
  country_df$pangea_id <- as.character(country_df$pangea_id)
  country_df$pangea_id <- as.factor(country_df$pangea_id)

  dna_seq <- country_df[c("pangea_id", "shiver_consensus")]
  dna_seq$shiver_consensus <- as.character(dna_seq$shiver_consensus)
  dna_seq$pangea_id <- as.character(dna_seq$pangea_id)

  # function requires to change name of columns
  colnames(dna_seq) <- c("seq.name", "seq.text")

  phylotools::dat2fasta(dna_seq, outfile = outfile_name)
}


#' Title Root phylogenetic tree and drop outgroup sequences
#'
#' Function to root a phylogenetic tree and drop outgroup sequences from the
#' tree.
#'
#' @param tree Phylogenetic tree of class phylo.
#' @param outgroups vector of string with names for outgroup sequences
#'
#' @return Rooted phylogenetic tree without outgroup sequences
#' @export
#'
#' @examples
#' # TO DO
root_and_drop_tips <- function(tree, outgroups){
  node <- ape::getMRCA(tree, outgroups)
  res <- try(phytools::reroot(tree, node.number = node, edgelabel = TRUE),
             silent = TRUE)

  if(class(res) == "try-error"){
    tree <- phytools::midpoint.root(tree)
    node <- ape::getMRCA(tree, outgroups)
    tree.r <- phytools::reroot(tree, node.number = node, edgelabel = TRUE)
    tree.r.new <- ape::drop.tip(tree.r, tip=outgroups)
  }else{
    tree.r <- phytools::reroot(tree, node.number = node, edgelabel = TRUE)
    tree.r.new <- ape::drop.tip(tree.r, tip=outgroups)
  }
  return(tree.r.new)
}



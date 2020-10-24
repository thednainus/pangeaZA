library(ape)
library(treedater)

#' Prunes redundant close global reference sequences for subsequent
#' coalescent analysis. 
#'
#' A CGR sequence is removed if it's closest 
#' neighbor in the tree is also a CGR. Thus the closest neighborr
#' for each CGR must be a non-CGR sequence
#' 
#' @param dtr A rooted phylogenetic tree of class ape::phylo or treedater. CGR lineages and *only* CGR sequencs must include the 'CGR' string in the tip label.
#' @value A pruned ape::phylo
#' @examples
#' load( system.file( 'data/env_treedater.rda', package = 'pangea') )
#' rrcgr( dtr_env2 ) -> tr 
#'
#' @export  
rrcgr <- function(dtr)
{
	tr <- dtr; class(tr) <- 'phylo' 
	iscgr <- grepl( tr$tip.label,  pattern = 'CGR' )
	desc <- lapply( 1:(Ntip(tr)+Nnode(tr)), function(i) c() )
	for ( i in 1:Ntip(tr))
		desc[[i]] <- i 
	allcgr <- rep( FALSE, Ntip(tr) + Nnode(tr))
	allcgr[which( iscgr )] <- TRUE
	
	ie <- postorder( tr )
	for ( i in ie){
		a <- tr$edge[i,1]
		u <- tr$edge[i,2] 
		desc[[a]] <- c( desc[[a]], desc[[u]] )
		allcgr[a] <- all ( iscgr[desc[[a]] ]) 
	}
	
	keepcgr <- c() 
	for ( i in ie){
		a <- tr$edge[i,1]
		u <- tr$edge[i,2] 
		if ( !allcgr[a] & allcgr[u] ){
			keepcgr <- c( keepcgr, desc[[u]][1] )
		}
	}
	dropcgr <- setdiff( tr$tip.label[iscgr], tr$tip.label[keepcgr] )
	tr2 <- drop.tip( tr, dropcgr )
	
	tr2
}

#' remove redundant cgr
#'
#' this version also removes cgrs with paraphyletic relationship 
#'
#' @export
rrcgr2 <- function(dtr)
{
	tr <- dtr; class(tr) <- 'phylo' 
	iscgr <- grepl( tr$tip.label,  pattern = 'CGR' )
	desc <- lapply( 1:(Ntip(tr)+Nnode(tr)), function(i) c() )
	for ( i in 1:Ntip(tr))
		desc[[i]] <- i 
	allcgr <- rep( FALSE, Ntip(tr) + Nnode(tr))
	allcgr[which( iscgr )] <- TRUE
	somecgr <- rep( FALSE, Ntip(tr) + Nnode(tr))
	somecgr[which(iscgr)] <- TRUE
	
	ie <- postorder( tr )
	for ( i in ie){
		a <- tr$edge[i,1]
		u <- tr$edge[i,2] 
		desc[[a]] <- c( desc[[a]], desc[[u]] )
		allcgr[a] <- all ( iscgr[desc[[a]] ]) 
		somecgr[a] <- any( iscgr[desc[[a]] ])
	}
	
	keepcgr <- c() 
	for( a in (Ntip(tr)+1):(Ntip(tr)+Nnode(tr))){
		uv <- tr$edge[ tr$edge[,1]==a ,2]
		u <- uv[1]
		v <- uv[2] 
		if ( !somecgr[u] & allcgr[v] ) {
			keepcgr <- c( keepcgr, desc[[v]][1] )
		}
		if ( allcgr[u] & !somecgr[v]){
			keepcgr <- c( keepcgr, desc[[u]][1] )
		}
	}
	dropcgr <- setdiff( tr$tip.label[iscgr], tr$tip.label[keepcgr] )
	tr2 <- drop.tip( tr, dropcgr )
	
	tr2
}

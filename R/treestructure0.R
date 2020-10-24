#' partition the tree into a set of clades within the desired size range 
#'
#' @export 
getsubtres <- function( tr , minsize = 200, maxsize = 500 , level = .25, NCPU = 1) 
{
	ts <- trestruct( tr, ncpu = NCPU , minCladeSize=200, level = level )
	tres <- lapply( ts$partitionSets, function(ps) keep.tip( tr, ps ))
	n <- sapply( tres, Ntip ) 
	list( 
	 justright = tres[ n <= maxsize ]
	 , toobig = tres [ n > maxsize ]
	)
}


#' partition the tree into a set of clades within the desired size range 
#'
#' If initial partitioning step produces clades that are too big, will 
#' continue partitioning at hiher sign levels until 
#' desired size is reached 
#'
#' @export 
getsubtresRecursively <- function( tr , minsize = 200, maxsize = 500 , level = .05, NCPU  = 1, incrementLevel = .15)
{
	done <- list( )
	todo <- list( tr ) 
	.level <- level - incrementLevel
	while( length( todo ) > 0 ){
		.level <- .level + incrementLevel  
		print( '------------------------------------' )
		print( c( .level, length( done ) , length(todo ) ) )
		
		ltodo <- length(todo )
		for (k in 1:ltodo){
			.tr <- todo[[k]]
			todo[[k]] <- FALSE
			o <- getsubtres( .tr, minsize, maxsize, level = .level   )
			done <- c( done, o$justright )
			todo <- c( todo , o$toobig )
		}
		todo <- todo [ sapply( todo, function(x) !is.logical(x) ) ]
		if ( .level > 1 )
			break
	}
	done 
}

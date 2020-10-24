

#' Generate tables and plots for fits of model0
#'
#' @param infns vector of filenames with MCMC fit in RDS format
#' @param burninpc proportion of mcmc to discard 
#' @param ncpu 
#' @param outdir where to store results; defaults to summary0
#' @param numSamples number of tractories to simulate
#' @value posterior sample (invisible)
#' @export 
summary_model0 <- function( infns , art_start, burninpc = .5,  ncpu = 1, outdir = NULL, numSamples = 100  ){
	stopifnot( require(phydynR) )
	library( BayesianTools )
	
	if (is.null( outdir )){
		outdir <- 'summary0' 
	}
	suppressWarnings( dir.create ( outdir ) )
	outdir <- paste0( outdir, '/' )
	
	chs <- createMcmcSamplerList( lapply( infns, readRDS ) )
	sink( file = paste0( outdir, 'mcmcSummary.txt' ))
	print( summary( chs ))
	sink() 
	
	O = getSample( chs )
	start <- floor( burninpc * nrow(O) )
	o <- O[ (start:nrow(O)), ]
#~ browser() 
	o <- o[ sample(1:nrow(o), size = numSamples, replace=FALSE), ]
	
	mod0 <- generate_model0(art_start=art_start)
	
	compute.traj0 <-  function(theta )
	{
		print(theta)
		p <- mod0$parms 
		p[ names(theta)] <- theta 
		x0 <- c(  src = p$src1980, i = p$i1980 , z =0) # NOTE this is set in compute.traj
		s <-  mod0$dm(p,  x0 , t0 = 1980, t1 =2015, res = 100 )
		s
	}
	i1980names <- colnames(O)[ grepl( 'i1980', colnames(O)) ]
	ntres <- length( i1980names )
	add.trajs <- function(ss){
		s <- ss[[1]]
		for (k in 2:length( ss )){
			.s <- ss[[k]] 
			s[[2]] <- lapply( 1:length( s[[2]]) , function(j) s[[2]][[j]] + .s[[2]][[j]] )
			s[[3]] <- lapply( 1:length( s[[3]]) , function(j) s[[3]][[j]] + .s[[3]][[j]] )
			s[[4]] <- lapply( 1:length( s[[3]]) , function(j) s[[4]][[j]] + .s[[4]][[j]] )
			s[[5]][, -1] <- s[[5]][, -1] + .s[[5]][, -1]
		}
		#browser() 
		s
	}
	compute.traj <- function( theta0){
		traj <- list() 
		for ( k in 1:ntres){
			theta <-  theta0[setdiff( names(theta0), i1980names ) ]
			theta <- c( theta , i1980 =  unname( theta0[ i1980names[k] ] ) )
			traj[[k]] <- compute.traj0( theta )
		}
		add.trajs( traj )
	}
	
	tfgys <- parallel::mclapply( 1:nrow(o), function(i){
		compute.traj( o[i, ] )
	}, mc.cores = ncpu)

	imat <- do.call( cbind, lapply( tfgys, function(s) s[[5]][, 'i' ] ) )
	zmat <- do.call( cbind, lapply( tfgys, function(s) s[[5]][, 'z' ] ) )
	qi <- t( sapply( 1:nrow(imat) , function(i) quantile( imat[i, ], prob = c( .5, .025, .975 )) ) )
	qz <- t( sapply( 1:nrow(imat) , function(i) quantile( zmat[i, ], prob = c( .5, .025, .975 )) ) )
	
	pzmat <- zmat / ( imat + zmat )
	qpz <- t( sapply( 1:nrow(pzmat) , function(i) quantile( pzmat[i, ], prob = c( .5, .025, .975 )) ) )

	time <- seq( 1980, 2015, length.out = 100 )
#~ browser()
	png( paste0( outdir, 'effinf.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qi, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'Effective infections')
	dev.off() 
	svg( paste0( outdir, 'effinf.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qi, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'Effective infections')
	dev.off() 
	
	png( paste0( outdir, 'ntreat.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'On ART')
	dev.off() 
	svg( paste0( outdir, 'ntreat.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'On ART')
	dev.off() 
	
	png( paste0( outdir, 'proptreat.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qpz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'Proportion on ART')
	dev.off() 
	svg( paste0( outdir, 'proptreat.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qpz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'Proportion on ART')
	dev.off() 
	
	newinf <- do.call( cbind, lapply( tfgys, function(s) {
		newinf <- sapply( s[[2]], function(FF) sum(FF[,2]) ) 
		rev( newinf  )
	}))
	qnewinf <- t( sapply( 1:nrow( newinf), function(i) quantile( newinf[i, ] , prob =c( .5, .025, .975 )))) 
	png( paste0( outdir, 'ninf.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qnewinf, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'New infections / year')
	dev.off() 
	svg( paste0( outdir, 'ninf.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qnewinf, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'New infections / year')
	dev.off() 
	
	
	propImport <- do.call( cbind, lapply( tfgys, function(s) {
		imp <- sapply( s[[3]], function(G) G[1,2] ) 
		crosstrans <-  sapply( s[[2]], function(FF) FF[1,2] )
		newinf <- sapply( s[[2]], function(FF) FF[2,2] ) 
		rev( (crosstrans + imp) / (crosstrans + imp + newinf ) )
	}))
	#~better as a box/whisker plot
	qpi <- t( sapply( 1:nrow( propImport), function(i) quantile( propImport[i, ] , prob =c( .5, .025, .975 )))) 
	sink( file = paste0( outdir, 'mcmcSummary.txt' ), append = TRUE )
	print('')
	print ('Proportion imports' )
	print( tail( qpi, 1)  )
	sink() 
	#X11(); matplot( time, qpi, type = 'l' , lty = c( 1, 3, 3), col = 'black' )
	
	bo <- lapply( mod0$betaNames, function(bn) o[, bn]  )
	png(paste0( outdir, 'beta.png')  , width = 4, height = 3 , units = 'in', res = 300 )
	par(mai = c( .5, .85, .3, .25) )
	boxplot( bo, names = mod0$betaTimes, ylab = 'Transmission rate' )
	dev.off() 
	svg(paste0( outdir, 'beta.svg')  , width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	boxplot( bo, names = mod0$betaTimes, ylab = 'Transmission rate' )
	dev.off() 
	
	invisible( o )
}










#' Make tfgys for multiple sets of input files  
#'
#' @export 
trajectories_model0 <- function( infns , art_start, mod0=NULL, burninpc = .5,  ncpu = 1, numSamples = 100 , onlyBestChain = FALSE ){
	stopifnot( require(phydynR) )
	library( BayesianTools )
	
	chss <- lapply( infns, readRDS )
	if (onlyBestChain & length(infns) > 1){
		 maps <- sapply( chss, function(ch)  MAP( ch )$valuesMAP[1] )
		 chss <- list( chss[[ which.max( maps ) ]] )
	}
	chs <- createMcmcSamplerList( chss )
	
	
#~ 	O = getSample( chs )
#~ 	O1 = getSample( chs, start = 15e3, thin = 10, parametersOnly=FALSE )
#~ 	start <- floor( burninpc * nrow(O) )
#~ 	o <- O[ (start:nrow(O)), ]
#~ 	o <- o[ sample(1:nrow(o), size = numSamples, replace=FALSE), ]
	
	O <- getSample( chs, numSamples = floor( numSamples / (1-burninpc) ) , parametersOnly = FALSE) 
	o <- tail( O , numSamples )
#~ browser() 
	
	if ( is.null(mod0)){
		warning('model unspecified; using default model0')
		mod0 <- generate_model0(art_start=art_start)
	}
	
	compute.traj0 <-  function(theta )
	{
		print(theta)
		p <- mod0$parms 
		p[ names(theta)] <- theta 
		x0 <- c(  src = p$src1980, i = p$i1980 , z =0) # NOTE this is set in compute.traj
		s <-  mod0$dm(p,  x0 , t0 = 1980, t1 =2015, res = 100 )
		s
	}
	i1980names <- colnames(O)[ grepl( 'i1980', colnames(O)) ]
	ntres <- length( i1980names )
	add.trajs <- function(ss){
		s <- ss[[1]]
		if ( length( ss) > 1 ){
			for (k in 2:length( ss )){
				.s <- ss[[k]] 
				s[[2]] <- lapply( 1:length( s[[2]]) , function(j) s[[2]][[j]] + .s[[2]][[j]] )
				s[[3]] <- lapply( 1:length( s[[3]]) , function(j) s[[3]][[j]] + .s[[3]][[j]] )
				s[[4]] <- lapply( 1:length( s[[3]]) , function(j) s[[4]][[j]] + .s[[4]][[j]] )
				s[[5]][, -1] <- s[[5]][, -1] + .s[[5]][, -1]
			}
		}
		#browser() 
		s
	}
	compute.traj <- function( theta0){
		traj <- list() 
		for ( k in 1:ntres){
			theta <-  theta0[setdiff( names(theta0), i1980names ) ]
			theta <- c( theta , i1980 =  unname( theta0[ i1980names[k] ] ) )
			traj[[k]] <- compute.traj0( theta )
		}
		add.trajs( traj )
	}
	
	tfgys <- parallel::mclapply( 1:nrow(o), function(i){
		compute.traj( o[i, ] )
	}, mc.cores = ncpu)
	
	list( tfgys = tfgys, o = o[, !grepl(pattern='i1980', colnames(o)) ])
}

#' @export 
trajectories_model0_multigene <- function( infns_list, ... )
{
	tfgyos <- lapply( infns_list, function(infn  ) 
	  trajectories_model0( infn , ... )
	 )
	tfgys <- do.call( c, lapply( tfgyos, '[[', 1 ) )
	o <- do.call( rbind, lapply( tfgyos, '[[', 2 ) )
	list( tfgys, o )
}



#' @export 
summary_model0_tfgyos <- function(tfgys, o , art_start,  outdir = NULL, axislog='', lquant = .025, uquant = .975 )
{
	library( coda ) 
	if (is.null( outdir )){
		outdir <- 'summary0_tfgys' 
	}
	suppressWarnings( dir.create ( outdir ) )
	outdir <- paste0( outdir, '/' )
	
	sink( file = paste0( outdir, 'mcmcSummary.txt' ))
	print( summary( coda::as.mcmc( o ))  )
	sink() 
	
	
	imat <- do.call( cbind, lapply( tfgys, function(s) s[[5]][, 'i' ] ) )
	zmat <- do.call( cbind, lapply( tfgys, function(s) s[[5]][, 'z' ] ) )
	qi <- t( sapply( 1:nrow(imat) , function(i) quantile( imat[i, ], prob = c( .5, lquant, uquant )) ) )
	qtotaliz <- t( sapply( 1:nrow(imat) , function(i) quantile( zmat[i, ] + imat[i, ] , prob = c( .5, lquant, uquant )) ) )
	qz <- t( sapply( 1:nrow(imat) , function(i) quantile( zmat[i, ], prob = c( .5, lquant, uquant )) ) )
	
	pzmat <- zmat / ( imat + zmat )
	qpz <- t( sapply( 1:nrow(pzmat) , function(i) quantile( pzmat[i, ], prob = c( .5, lquant, uquant )) ) )

	#time <- seq( 1980, 2015, length.out = 35*6)
	time <- rev( tfgys[[1]][[1]] ) 
	
	
	# R(t) 
	M = mod0 <- generate_model0(art_start = art_start )
#~ 	Rt <- sapply( 1:nrow(o), function(i){
#~ 		p <- as.list( o[i,] )
#~ 		p$art_start = unname( art_start ) 
#~ 		b <- M$parms$beta.t( time, p )
#~ 		r <- sapply( time, function(tt) M$parms$rho.t( tt, p ) )
#~ 		b / ( r + M$parms$gamma + M$parms$mu )
#~ 	})
	Rt <- sapply( tfgys, m0R.t )
	qRt <- t( sapply( 1:nrow(Rt) , function(i) quantile( Rt[i, ], prob = c( .5, lquant, uquant )) ) )
	
#~ browser()
	png( paste0( outdir, 'R.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qRt, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'R(t)')
	abline( h = 1,col = 'red')
	dev.off() 
	svg( paste0( outdir, 'R.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qRt, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'R(t)')
	abline( h = 1,col = 'red')
	dev.off() 
	
	
	png( paste0( outdir, 'i.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qi, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'I(t)', log = axislog)
	dev.off() 
	svg( paste0( outdir, 'i.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qi, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'I(t)' , log = axislog)
	dev.off() 
	
		
		png( paste0( outdir, 'iztotal.png'), width = 4, height = 3 , units = 'in', res  = 300)
		par(mai = c( .5, .85, .3, .25) )
		matplot( time, qtotaliz, type = 'l', lty = c(1, 3, 3), col = 'black' 
		  , xlab = '', ylab = 'I(t) +  Z(t)', log = axislog)
		dev.off() 
		svg( paste0( outdir, 'iztotal.svg'), width = 4, height = 3 )
		par(mai = c( .5, .85, .3, .25) )
		matplot( time, qtotaliz, type = 'l', lty = c(1, 3, 3), col = 'black' 
		  , xlab = '', ylab = 'I(t) +  Z(t)', log = axislog)
		dev.off() 
	
	png( paste0( outdir, 'ntreat.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'On ART')
	dev.off() 
	svg( paste0( outdir, 'ntreat.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'On ART')
	dev.off() 
	
	png( paste0( outdir, 'proptreat.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qpz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'Proportion on ART')
	dev.off() 
	svg( paste0( outdir, 'proptreat.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qpz, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'Proportion on ART')
	dev.off() 
	
	newinf <- do.call( cbind, lapply( tfgys, function(s) {
		newinf <- sapply( s[[2]], function(FF) sum( FF[,2]) )
		rev( newinf  )
	}))
	qnewinf <- t( sapply( 1:nrow( newinf), function(i) quantile( newinf[i, ] , prob =c( .5, lquant, uquant )))) 
	png( paste0( outdir, 'ninf.png'), width = 4, height = 3 , units = 'in', res  = 300)
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qnewinf, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'New infections / year', log = axislog)
	dev.off() 
	svg( paste0( outdir, 'ninf.svg'), width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	matplot( time, qnewinf, type = 'l', lty = c(1, 3, 3), col = 'black' 
	  , xlab = '', ylab = 'New infections / year', log = axislog)
	dev.off() 
	
	
	propImport <- do.call( cbind, lapply( tfgys, function(s) {
		imp <- sapply( s[[3]], function(G) G[1,2] ) 
		crosstrans <- sapply( s[[2]], function(FF) FF[1,2] ) 
		newinf <- sapply( s[[2]], function(FF) FF[2,2] )
		rev( (crosstrans + imp) / (crosstrans + imp + newinf ) )
	}))
	#~better as a box/whisker plot
	qpi <- t( sapply( 1:nrow( propImport), function(i) quantile( propImport[i, ] , prob =c( .5, lquant, uquant )))) 
	sink( file = paste0( outdir, 'mcmcSummary.txt' ), append = TRUE )
	print('')
	print ('Proportion imports' )
	print( tail( qpi, 1)  )
	sink() 
	#X11(); matplot( time, qpi, type = 'l' , lty = c( 1, 3, 3), col = 'black' )
	
	bo <- lapply( mod0$betaNames, function(bn) o[, bn]  )
	png(paste0( outdir, 'beta.png')  , width = 4, height = 3 , units = 'in', res = 300 )
	par(mai = c( .5, .85, .3, .25) )
	boxplot( bo, names = mod0$betaTimes, ylab = 'Transmission rate' )
	dev.off() 
	svg(paste0( outdir, 'beta.svg')  , width = 4, height = 3 )
	par(mai = c( .5, .85, .3, .25) )
	boxplot( bo, names = mod0$betaTimes, ylab = 'Transmission rate' )
	dev.off() 
	
#~ 	browser()
	
	invisible( o )



}

#' Compute R(t) by integrating over future hazard of removal and transmissions 
#' 
#' @export 
m0R.t <- function ( s) {
	i <- 2
	y <- rev( sapply( s$sizes, '[', i ) )
	d <- rev( sapply( s$deaths, '[', i ))
	b <- rev( sapply( s$births, function(x) x[i,i] ))
	x <- rev( s$times )
	dx <- diff(x)[1] 
	
	r <- dx * b / y 
	r <- c( r, rep( tail(r, 1), 100 ))
	loghaz <- log (1- d*dx / y )
	lasthazard <- log (1- tail(d,1)*dx / tail(y,1) )
	loghaz <- c( loghaz, rep(lasthazard, 100 ) )
	logsurv <- cumsum( loghaz)
	R <- sapply( 1:length( x ), function(k){
		j <-  k:length(r)
		sum( exp( logsurv[j] - logsurv[k]  ) * r[j]  )
	})
	R
}


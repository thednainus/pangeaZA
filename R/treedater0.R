#' Does a treedater analysis for gag, env or pol and usign the south africa data 
#' 
#' @export 
sa_treedater <- function( TREEFN = 'sa_ml_gag.nwk'
	, SRCMDFN = 'SA_src_gag_metadata_mod.csv'  #~ SA_src_pol_metadata_mod.csv
	, PANGEAFN = "../../PANGEA_Extract_Volz_2019-04-16.csv"
	, OG = c("K03454_gag", "AY371157_gag")
	, seqlen = 1476 # the length of the HIV sequences
	, NCPU = 8
	, outstem = 'gag'
){
	
# Script to get data for gag gene to be used with treedater
library(lubridate)
#~ library(senegalHIVmodel)
library(ape)
library(pangea)

# Get data for South African sequences

#~ source("get_SA_data.R") #
	{
		
	pangea_data = pd <- read.csv(PANGEAFN ,
							na.strings=c("","NA"))
		
		pd$pangea_id <- as.character(pd$pangea_id)
		pd$pangea_id <- as.factor(pd$pangea_id)

		pd$geo_country <- as.character(pd$geo_country)
		pd$geo_country <- as.factor(pd$geo_country)

		pd$gender <- as.character(pd$gender)
		pd$gender <- as.factor(pd$gender)

		pd$ymd = ymd(pd$sample_date)

		pd$decimal <- decimal_date( pd$ymd )

		rownames( pd ) <- as.character( pd$pangea_id  )

		# get data from South Africa sequences
		sa <- subset(pd, main_cohort_id == "South Africa" | geo_country == "SouthAfrica")
		sa <- subset(sa, gender != "Unknown")

		# get pangea data that was used as CGR sequences
		load( system.file( "data/SA/SA_pangea_cgr.rda", package = 'pangea'))

		# get sample_date for each South African sequence
		sa_dates <- sa[c("pangea_id", "geo_country", "sample_date", "shiver_consensus", 'ymd' , 'decimal')]
		rownames( sa_dates ) <- as.character( sa_dates$pangea_id )

	}

# read rooted plylogenetic tree for gag
# Here the tree is unrooted and I will root it using R before droping
# outgroup sequences
tr <- read.tree( TREEFN  )

# Root phylogenetic tree and drop outgroup sequences
tr2 <- root_and_drop_tips(tree = tr,
                              outgroups = OG)


# get sequences that is in the phylogenetic tree for pangea CGR
# the lapply function will simply remove the "_CGR" from some phylogenetic tip
# label to match sequences in the pangea dataset.
pangea_gag <- pangea_data[pangea_data$pangea_id %in%
                            lapply(tr2$tip.label,
                                   function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1]),]
                                   
sts <- setNames( rep( NA, Ntip( tr2 )), tr2$tip.label )
x <- names( sts  ) [ names(sts) %in% rownames(sa_dates )]
sts[x] <- sa_dates[x, ]$decimal

## fill in pangea cgr's
x <- tr2$tip.label[ grepl( tr2$tip.label, pattern= 'CGR' ) ]
x <- sapply( strsplit( x, '_' ), '[', 1 )
x <- x[ x %in% pd$pangea_id ]
rownames( pd ) <- as.character( pd$pangea_id  )
z <- paste0( x, '_CGR' )
sts[z] <- pd[x,]$decimal



## missing 
mst <- data.frame( lower = 1980, upper = 2015, taxon = names(sts)[is.na(sts)] )
rownames( mst ) <- mst$taxon
mst <- mst[ , c('lower', 'upper') ]

## fill in genbank cgr's bounds 
#~ y <- tr2$tip.label[ grepl( tr2$tip.label, pattern= 'CGR' ) ]
#~ y <- sapply( strsplit( y, '_' ), '[', 1 )
#~ y <- setdiff( y, x )
#~ y <- sapply( strsplit(y , '\\.' ), '[', 1 )
y <- sapply( strsplit(rownames(mst) , '\\.' ), '[', 1 )
cgrdf <- read.csv(SRCMDFN )
cgrdf$sid <- as.character( cgrdf$GenBank )
for ( i in 1:nrow(mst)){
	yy <- y[i]
	z <- rownames( mst )[i]
	if ( yy %in% cgrdf$sid ){
		w <- as.character( cgrdf$Collection_year[ match( yy, cgrdf$sid ) ] )
		if ( nchar(w)==4){
			alb <- as.numeric(w)
			aub <- alb + 1
		} else if ( grepl( w, pattern = 'to' ) ){
			a <- strsplit( w, ' to ' )[[1]]
			if ( nchar( a[1] )==4){
				alb <- as.numeric( a[1] )
			} else{
				alb <- decimal_date( dmy( a[1] ) )
			}
			if ( nchar( a[2] )==4){
				aub <- as.numeric( a[2] )
			} else{
				aub <- decimal_date( dmy( a[2] ) )
			}
		} else {
			alb  <- decimal_date( dmy( w ) )
			aub <- alb  + 1/12
		}
		if ( !is.na( alb ) & !is.na( aub )){
			if ((aub < 2015) & (alb > 1980) & (alb < 2015) & (aub > 1980) & (alb <= aub)){
				mst[z,]$lower = alb
				mst[z,]$upper = aub
			}
		}
	}
}








(dtr <- dater(tr2, sts, seqlen, estimateSampleTimes = mst, ncpu = NCPU , quiet = FALSE))
ot0 <- outlierTips( dtr )
saveRDS(ot0, paste0( outstem, "_outliers.rds") )


tr3 <- drop.tip( tr2, as.character( ot0$taxon[ ot0$q < .05] ) )
dtr2 <- dater(tr3, sts, seqlen, estimateSampleTimes = mst, omega0 = dtr$meanRate , ncpu = NCPU )


saveRDS(dtr2, file = paste0(outstem, "_treedater.rds") )

dtr2
}

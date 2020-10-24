library(phydynR)

# median survival time, africa, age 30 at seroconversion (mangal et al PMCID: PMC5414573)
.SURVTIME <- 10.2
# life expectancy is 62.7 yrs
# median age of seroconversion is 30 years (mangal)
.NATMORT <- 1/(62.7-30) # nat mort


#' Generate a simple HIV epidemic model for use by phydynR
#'
#" This model has the following demes:
#' * src : global reservoir
#' * i : infected and infectious
#' * z : infected and on ART
#' Transmission occurs at a rate beta(t) i(t) where beta is a linear
#' interpolating spline. The treatment rate i->z is initially zero but
#' increases linearly to 2015
#'
#' @value A list with elements
#' * dm : a demographic process model
#' * parms: a default list of parameters
#' * x0 : a default vector of start conditions
#'
#' @export
generate_model0<- function(art_start = 2005)
{

parms <- list(
 srcGrowth = 0.035
 , src1980  =  1e4
 , i1980 = 1
 , gamma  = 1/.SURVTIME # disease mort
 , mu = .NATMORT
 , importRate = 1/20
 , rho2017 = 1/4.5 #  rate art initiation
 , beta1980 = 10 / .SURVTIME
 , beta1990 = 6 / .SURVTIME
 , beta2000 = 4 / .SURVTIME
 , beta2005 = 4 / .SURVTIME
 , beta2010 = 2 / .SURVTIME
 , beta2015 = 2 / .SURVTIME
 , art_start = art_start
)



betaTimes <- c( 1980, 1990, 2000, 2005, 2010, 2015 )
betaNames <- paste0( 'beta', betaTimes )

parms$beta.t <- function(t, p){
	approx( betaTimes, unlist( p[ betaNames ] ) , xout = t, rule = 2)$y
}

parms$rho.t <- function(t, p){
	if ( t < p$art_start )
	  return (0 )
	approx( c( p$art_start, 2017), c(0, p$rho2017), xout = t , rule = 2)$y
}

demes <- c( 'src', 'i' )
nondemes <- c( 'z' )

x0 <- c( src = 1e3, i = 1 , z = 0 )

eqns <- setup.model.equations( demes, nondemes ) #  cpy from senegalHIV
attach( eqns) #birth, deaths, migs

births['src','src'] <- ' src * parms$srcGrowth'
births['i', 'i' ] <- 'i * parms$beta.t(t, parms)'

migs[ 'src', 'i' ] <- 'i * parms$importRate'
migs['i', 'src'] <- 'i * parms$importRate'

deaths['i'] <- 'i * parms$mu  + i * parms$gamma + i * parms$rho.t( t, parms ) '

nonDemeDynamics['z'] <- ' i * parms$rho.t( t, parms ) - z * parms$mu '


dm <- build.demographic.process( births = births
 , migrations = migs
 , deaths = deaths
 , nonDemeDynamics = nonDemeDynamics
 , parameterNames = names(parms )
 , sde = FALSE
 , rcpp = FALSE
)

list( dm = dm
, parms = parms
, x0 = x0
, betaTimes = betaTimes
, betaNames = betaNames
, m = m
, mm = mm
, demes = demes
, nondemes = nondemes
)

}

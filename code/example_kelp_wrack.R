## This script implements a worked empirical example based on a subset of the 
## results in Walter et al. (2024) Spatial synchrony cascades across ecosystem 
## boundaries and up food webs via resource subsidies. PNAS 121(2) e23152120.

## The data provided to accompany this example have been pre-processed to 
## (where appropriate) spatially aggregate, organize, and gap-fill.
## Original data sources are publicly accessible:

## Santa Barbara Coastal LTER, J. Dugan, SBC LTER: Beach: 
## Time-series of beach wrack cover and biomass, ongoing since 2008 ver 16. 
## Environmental Data Initiative (2021). 
## https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.40.16.

## T. Bell, SBC LTER: Reef: California kelp canopy and environmental variable dynamics ver 1. 
## Environmental Data Initiative (2023). https://doi.org/10.6073/pasta/c40db2c8629cfa3fbe80fdc9e086a9aa. 


## Examine coherence between kelp forest production and kelp wrack deposition

rm(list=ls())
graphics.off()

library(wsyn)


#This sets the working directory to the ./code subdirectory but only works if using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

## load datasets ----------------------------------------------------------------------------------

#these files have already been organized into locations-by-time matrices and been gap-filled by
#replacing NAs with month-specific medians
wrack <- as.matrix(read.csv("../data/wrack.csv"))
beachwidth <- as.matrix(read.csv("../data/beachwidth.csv"))
wavediff <- as.matrix(read.csv("../data/wavediff.csv"))
kelp <- as.matrix(read.csv("../data/kelpproduction.csv"))

#replicate replicate the region-wide kelp time series 5x to match site dimensions
kelp <- rbind(c(kelp),c(kelp),c(kelp),c(kelp),c(kelp))

#this is the time (year and month) information corresponding to columns in the data matrices
yearmonth <- read.csv("../data/time.csv")
tt=1:132 #vector of timesteps


## apply cleandat() to normalize, detrend, and scale ----------------------------------------------

beachwidth.cln <- cleandat(beachwidth, tt, clev=5)$cdat
wavediff.cln <- cleandat(wavediff, tt, clev=5)$cdat
kelp.cln <- cleandat(kelp, tt, clev=5)$cdat
wrack.cln <- cleandat(wrack, tt, clev=5)$cdat


## Create wavelet mean fields ---------------------------------------------------------------------


wmf.wrack <- wmf(wrack.cln, tt)
wpmf.wrack <- wpmf(wrack.cln, tt, sigmethod="quick")

plotmag(wmf.wrack)
plotmag(wpmf.wrack)


## Wavelet spatial coherence analyses -------------------------------------------------------------

#set timescale ranges
subann <- c(0,8) #months
annual <- c(8,16) #months
intann <- c(16,60) #months

#nrand is relatively small to speed computation. Often 10,000 is used for publications with the fast method
wrackXkelp <- coh(wrack.cln, kelp.cln, tt, sigmethod="fast", norm="powall", nrand=1000)
wrackXkelp <- bandtest(wrackXkelp, subann)
wrackXkelp <- bandtest(wrackXkelp, annual)
wrackXkelp <- bandtest(wrackXkelp, intann)
plotmag(wrackXkelp) #coherence on subannual, annual and interannual 
plotphase(wrackXkelp) #phase lagged on annual and interannial 

wrackXwidth <- coh(wrack.cln, beachwidth.cln, tt, sigmethod="fast", norm="powall", nrand=1000)
wrackXwidth <- bandtest(wrackXwidth, subann)
wrackXwidth <- bandtest(wrackXwidth, annual)
wrackXwidth <- bandtest(wrackXwidth, intann)
plotmag(wrackXwidth) #coherence on seasonal, annual, interannual
plotphase(wrackXwidth) #in phase on annual, borderline in-phase or lagged on long timescales

wrackXwavediff <- coh(wrack.cln, wavediff.cln, tt, sigmethod="fast", norm="powall", nrand=1000)
wrackXwavediff <- bandtest(wrackXwavediff, subann)
wrackXwavediff <- bandtest(wrackXwavediff, annual)
wrackXwavediff <- bandtest(wrackXwavediff, intann)
plotmag(wrackXwavediff) #coherence on seasonal, annual, interannual
plotphase(wrackXwavediff) #in phase on annual, borderline in-phase or lagged on long timescales



## Wavelet linear model of wrack ------------------------------------------------------------------

#makes a wavelet linear model based on the set of putative predictors that were coherent with
#wrack on the annual timescale band

datlist <- list(wrack.cln, kelp.cln, wavediff.cln, beachwidth.cln)
wlm.annual <- wlm(datlist, tt, resp=1, pred=c(2,3,4), norm="powall", f0=1)
#get whole-model p-value by 'dropping' all predictors
wlmtest.annual <- wlmtest(wlm.annual, drop=2:4, sigmethod = "fft", nrand=100)
wlmtest.annual <- bandtest(wlmtest.annual, annual)
print(wlmtest.annual$bandp)

#note: can drop a single predictor to get a p-value for that predictor

#look at synchrony explained across all timescales
syncexpl.annual <- syncexpl(wlm.annual)
print(syncexpl.annual)

#subset to wavelet components in annual timescale band
syncexpl.annual <- syncexpl.annual[syncexpl.annual$timescales > min(annual) 
                                   & syncexpl.annual$timescales < max(annual),]
print(syncexpl.annual)

#average across whole annual timescale band
print(colMeans(syncexpl.annual))

#overall synchrony explained in the annual timescale band
print(mean(syncexpl.annual$syncexpl)/mean(syncexpl.annual$sync))

#overall 'crossterms' across the annual timescale band
mean(syncexpl.annual$crossterms)/mean(syncexpl.annual$syncexpl)

#crossterms are a measure of spatial independence that evaluates 
#an assumption of the synchrony explained partitioning scheme.
#We have considered values <10% of the synchrony explained acceptable.

#output and plot predicted synchrony

predsync.annual <- predsync(wlm.annual)


plotmag(predsync.annual, zlims=c(0,max(Mod(wmf.wrack$values), na.rm=TRUE)))
abline(h=log2(min(annual)), lwd=2)
abline(h=log2(max(annual)), lwd=2)


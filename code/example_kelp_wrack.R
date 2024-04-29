### EDITING IN PROGRESS



## Examine coherence between kelp forest production and kelp wrack deposition

rm(list=ls())

library(wsyn)
#library(viridis)
# library(rgdal)
# library(ncf)
# library(raster)
# library(lubridate)
# library(car)
# library(fields)
# library(RColorBrewer)

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
wavediff.cln <- cleandat(wavesdiff, tt, clev=5)$cdat
kelp.cln <- cleandat(kelp, tt, clev=5)$cdat
wrack.cln <- cleandat(wrack, tt, clev=5)$cdat


#set timescale ranges
subann <- c(0,8) #months
annual <- c(8,16) #months
intann <- c(16,60) #months





## CUSTOM PLOTTING FUNCTION -----------------------------------------------------------------------

plotwmf<-function(object,zlims=NULL,neat=TRUE,colorfill=NULL,colorbar=TRUE,title=NULL,filename=NA,xlocs=NULL,xlabs=NULL,...)
{
  wav<-Mod(get_values(object))
  times<-get_times(object)
  timescales<-get_timescales(object)
  
  if(is.null(zlims)){
    zlims<-range(wav,na.rm=T)
  }else
  {
    rg<-range(wav,na.rm=T)
    if (rg[1]<zlims[1] || rg[2]>zlims[2])
    {
      stop("Error in plotmag.tts: zlims must encompass the z axis range of what is being plotted")
    }
  }
  if(neat){
    inds<-which(!is.na(colMeans(wav,na.rm=T)))
    wav<-wav[,inds]
    timescales<-timescales[inds]
  }
  if(is.null(colorfill)){
    jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    colorfill<-grDevices::colorRampPalette(jetcolors)
  }
  ylocs<-pretty(timescales,n=8)
  if(is.null(xlocs)){
    xlocs<-pretty(times,n=8)
  }
  
  if(is.null(xlabs)){
    xlabs<-xlocs
  }
  
  if (!is.na(filename))
  {
    grDevices::pdf(paste0(filename,".pdf"))
  }
  if (!colorbar)
  {
    graphics::image(x=times,y=log2(timescales),z=wav,xlab="Time",zlim=zlims,
                    ylab="Timescale",axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlabs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }else
  {
    fields::image.plot(x=times,y=log2(timescales),z=wav,xlab="Time",zlim=zlims,
                       ylab="Timescale",axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlabs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }
  if (!is.na(filename))
  {
    grDevices::dev.off()
  }
}


## WAVELET MEAN FIELDS ----------------------------------------------------------------------------

wmf.nearprod <- wmf(near.production.cln, tt, f0=1, scale.max.input = 60)
wmf.wrack <- wmf(wrack.bysite.cln, tt, f0=1, scale.max.input = 60)

wpmf.wrack <- wpmf(wrack.bysite.cln, tt, sigmethod="quick", f0=1, scale.max.input = 60)
wpmf.nearprod <- wpmf(near.production.cln, tt, sigmethod = "quick", f0=1, scale.max.input = 60)
wpmf.width <- wpmf(width.cln, tt, sigmethod="quick", f0=1, scale.max.input = 60)



## WAVELET SPATIAL COHERENCE: WRACK DRIVERS -------------------------------------------------------

wrackXnearProd <- coh(wrack.bysite.cln, near.production.cln, tt, sigmethod="fast", norm="powall", nrand=10000, f0=1)
wrackXnearProd <- bandtest(wrackXnearProd, subann)
wrackXnearProd <- bandtest(wrackXnearProd, annual)
wrackXnearProd <- bandtest(wrackXnearProd, intann)
plotmag(wrackXnearProd) #no coherence 
plotphase(wrackXnearProd) #antiphase on seasonal timescales

wrackXallProd <- coh(wrack.bysite.cln, all.production.cln, tt, sigmethod="fast", norm="powall", nrand=10000, f0=1)
wrackXallProd <- bandtest(wrackXallProd, subann)
wrackXallProd <- bandtest(wrackXallProd, annual)
wrackXallProd <- bandtest(wrackXallProd, intann)
plotmag(wrackXallProd) #coherence on subannual, annual and interannual 
plotphase(wrackXallProd) #phase lagged on annual and interannial 

wrackXwidth <- coh(wrack.bysite.cln, width.cln, tt, sigmethod="fast", norm="powall", nrand=10000, f0=1)
wrackXwidth <- bandtest(wrackXwidth, subann)
wrackXwidth <- bandtest(wrackXwidth, annual)
wrackXwidth <- bandtest(wrackXwidth, intann)
plotmag(wrackXwidth) #coherence on seasonal, annual, interannual
plotphase(wrackXwidth) #in phase on annual, borderline in-phase or lagged on long timescales

wrackXwavediff <- coh(wrack.bysite.cln, wavediff.cln, tt, sigmethod="fast", norm="powall", nrand=10000, f0=1)
wrackXwavediff <- bandtest(wrackXwavediff, subann)
wrackXwavediff <- bandtest(wrackXwavediff, annual)
wrackXwavediff <- bandtest(wrackXwavediff, intann)
plotmag(wrackXwavediff) #coherence on seasonal, annual, interannual
plotphase(wrackXwavediff) #in phase on annual, borderline in-phase or lagged on long timescales



## WAVELET LINEAR MODELS: WRACK DRIVERS -----------------------------------------------------------

timescales <- wpmf.wrack$timescales
datlist_wrack <- list(wrack.bysite.cln, near.production.cln, all.production.cln, wavediff.cln, width.cln)

wlm.subann <- wlm(datlist_wrack, tt, resp=1, pred=c(3,4,5), norm="powall", f0=1)
wlmtest.subann <- wlmtest(wlm.subann, drop=2:3, sigmethod = "fft", nrand=10000)
wlmtest.subann <- bandtest(wlmtest.subann, subann)
wlmtest.subann$bandp

syncexpl.subann <- syncexpl(wlm.subann)
colMeans(syncexpl.subann[timescales > min(subann) & timescales < max(subann),])
Arg(colMeans(wlm.subann$coefs[timescales > min(subann) & timescales < max(subann),]))/pi

mean(syncexpl.subann$syncexpl[timescales > min(subann) & timescales < max(subann)])/
  mean(syncexpl.subann$sync[timescales > min(subann) & timescales < max(subann)])
mean(syncexpl.subann$crossterms[timescales > min(subann) & timescales < max(subann)])/
  mean(syncexpl.subann$syncexpl[timescales > min(subann) & timescales < max(subann)])


wlm.annual <- wlm(datlist_wrack, tt, resp=1, pred=c(3,4,5), norm="powall", f0=1)
wlmtest.annual <- wlmtest(wlm.annual, drop=2:4, sigmethod = "fft", nrand=10000)
wlmtest.annual <- bandtest(wlmtest.annual, annual)
wlmtest.annual$bandp

syncexpl.annual <- syncexpl(wlm.annual)
colMeans(syncexpl.annual[timescales > min(annual) & timescales < max(annual),])
Arg(colMeans(wlm.annual$coefs[timescales > min(annual) & timescales < max(annual),]))/pi

mean(syncexpl.annual$syncexpl[timescales > min(annual) & timescales < max(annual)])/
  mean(syncexpl.annual$sync[timescales > min(annual) & timescales < max(annual)])
mean(syncexpl.annual$crossterms[timescales > min(annual) & timescales < max(annual)])/
  mean(syncexpl.annual$syncexpl[timescales > min(annual) & timescales < max(annual)])


wlm.intann <- wlm(datlist_wrack, tt, resp=1, pred=c(3,5), norm="powall", f0=1)
wlmtest.intann <- wlmtest(wlm.intann, drop=2:3, sigmethod = "fft", nrand=10000)
wlmtest.intann <- bandtest(wlmtest.intann, intann)
wlmtest.intann$bandp

bandtest(wlmtest.intann, c(30,60))$bandp

syncexpl.intann <- syncexpl(wlm.intann)
colMeans(syncexpl.subann[timescales > min(annual) & timescales < max(annual),])
Arg(colMeans(wlm.intann$coefs[timescales > min(intann) & timescales < max(intann),], na.rm=T))/pi

mean(syncexpl.intann$syncexpl[timescales > min(intann) & timescales < max(intann)], na.rm=T)/
  mean(syncexpl.intann$sync[timescales > min(intann) & timescales < max(intann)], na.rm=T)
mean(syncexpl.intann$crossterms[timescales > min(intann) & timescales < max(intann)], na.rm=T)/
  mean(syncexpl.intann$syncexpl[timescales > min(intann) & timescales < max(intann)], na.rm=T)



## Figure demonstrating synchrony explained by drivers --------------------------------------------


predsync.subann <- predsync(wlm(datlist_wrack,tt, resp=1, pred=c(3,4,5), norm="powall", f0=1)) #waves, beach width
predsync.annual <- predsync(wlm(datlist_wrack,tt, resp=1, pred=c(3,4,5), norm="powall", f0=1)) #all production, waves, width
predsync.intann <- predsync(wlm(datlist_wrack,tt, resp=1, pred=c(3,5), norm="powall", f0=1)) #all production, width

predsync_comb <- predsync.subann
predsync_comb$values[!is.infinite(predsync_comb$values)] <- NA #clear values
predsync_comb$values[,predsync_comb$timescales>=min(subann) & predsync_comb$timescales<max(subann)] <- predsync.subann$values[,predsync_comb$timescales>=min(subann) & predsync_comb$timescales<max(subann)]
predsync_comb$values[,predsync_comb$timescales>min(annual) & predsync_comb$timescales<max(annual)] <- predsync.annual$values[,predsync_comb$timescales>min(annual) & predsync_comb$timescales<max(annual)]
predsync_comb$values[,predsync_comb$timescales>min(intann) & predsync_comb$timescales<max(intann)] <- predsync.intann$values[,predsync_comb$timescales>min(intann) & predsync_comb$timescales<max(intann)]


laymat <- matrix(1:8, nrow=4, ncol=2, byrow=FALSE)
xlabs <- c("2010-1","2012-1","2014-1","2016-1","2018-1")
jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
colorfill<-grDevices::colorRampPalette(jetcolors)
pal <- brewer.pal(5,"Set3")


altsqrt <- function(x){
  out <- rep(NA, length(x))
  out[x >= 0] <- sqrt(x[x >= 0])
  out[x < 0] <- -1*sqrt(abs(x[x < 0]))
  return(out)
}

ts1<-syncexpl.subann$timescales>=min(subann) & syncexpl.subann$timescales<max(subann)
ts2<-syncexpl.subann$timescales>=min(annual) & syncexpl.subann$timescales<max(annual)
ts3<-syncexpl.subann$timescales>=min(intann) & syncexpl.subann$timescales<max(intann)

y1 <- colMeans(syncexpl.subann[ts1,c("pred1","pred2","pred3","interactions","syncexpl")]/syncexpl.subann$sync[ts1])*100
y2 <- colMeans(syncexpl.annual[ts2,c("pred1","pred2","pred3","interactions","syncexpl")]/syncexpl.annual$sync[ts2])*100
y3 <- colMeans(syncexpl.intann[ts3,c("pred1","pred2","interactions","syncexpl")]/syncexpl.intann$sync[ts3])*100

wrack.syncexpl.bandavg <- cbind(c(y1)
                                ,y2
                                ,c(y3[1],NA,y3[2:4])
)


laymat <- matrix(1:10, nrow=5, byrow=FALSE)
pal <- brewer.pal(5,"Set3")


pdf("./figures/fig_wrackSynchrony_ann8_16_alt.pdf", width=3.5, height=7)

layout(laymat, widths=c(0.9,0.1), heights=c(1,1,0.85,0.85,1.15))

par(mar=c(2.1,1.8,1.6,0.1), oma=c(0, 2.5, 0, 0.4), mgp=c(2,0.8,0))

plotwmf(wmf.wrack, xlocs=seq(from=12, by=24, length.out=5), xlabs=xlabs, colorbar=FALSE)
q<-stats::quantile(wpmf.wrack$signif[[2]], 0.95)
graphics::contour(x=wpmf.wrack$times,y=log2(wpmf.wrack$timescales),z=Mod(wpmf.wrack$values),levels=q,drawlabels=F,lwd=1,
                  xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F, add=T)
mtext("Wrack synchrony", cex=0.8, line=0.2)
abline(h=log2(max(subann)), col="white", lwd=1.5)
abline(h=log2(max(annual)), col="white", lwd=1.5)
mtext("a)", cex=0.8, at=0, line=0.2)

plotwmf(predsync_comb, xlocs=seq(from=12, by=24, length.out=5), zlim=c(0,1.6), xlabs=xlabs, neat=FALSE, colorbar=FALSE)
mtext("Predicted synchrony", cex=0.8, line=0.2)
abline(h=log2(max(subann)), col="white", lwd=1.5)
abline(h=log2(max(annual)), col="white", lwd=1.5)
mtext("b)", cex=0.8, at=0, line=0.2)

par(mar=c(2.8,1.8,1.6,0.1), mgp=c(1.8,0.8,0))
plot(NA, NA, ylim=c(0,1.2), xlim=range(tt), ylab="Synchrony", xlab="", xaxt="n", xaxs="i")
lines(tt, Mod(wmf.wrack$values[,which.min(abs(timescales-12))]))
lines(tt, Mod(predsync_comb$values[,which.min(abs(timescales-12))]), col=pal[5])
axis(1, at=seq(from=12, by=24, length.out=5), labels=xlabs)
legend("bottomright", bty="n", lty=1, col=c("black",pal[5]), legend=c("Synchrony","Pred. sync."))
mtext("c)", cex=0.8, at=1, line=0.2)

par(mar=c(2.8,1.8,1.6,0.1), mgp=c(1.8,0.8,0))
plot(NA, NA, ylim=c(0,0.8), xlim=log2(c(2,60)), ylab="Synchrony", xlab="Timescale (months)", xaxt="n")
axis(1, at=log2(c(2,4,8,16,24,36,48,60)), labels=c(2,4,8,16,24,36,48,60))
lines(log2(syncexpl.subann$timescales), syncexpl.subann$sync)
lines(log2(syncexpl.subann$timescales[syncexpl.subann$timescales>=min(subann) & syncexpl.subann$timescales<max(subann)]),
      syncexpl.subann$syncexpl[syncexpl.subann$timescales>=min(subann) & syncexpl.subann$timescales<max(subann)], col=pal[5])
lines(log2(syncexpl.annual$timescales[syncexpl.annual$timescales>min(annual) & syncexpl.annual$timescales<max(annual)]),
      syncexpl.annual$syncexpl[syncexpl.annual$timescales>min(annual) & syncexpl.annual$timescales<max(annual)], col=pal[5])
lines(log2(syncexpl.intann$timescales[syncexpl.intann$timescales>min(intann) & syncexpl.intann$timescales<max(intann)+1]),
      syncexpl.intann$syncexpl[syncexpl.intann$timescales>min(intann) & syncexpl.intann$timescales<max(intann)+1], col=pal[5])
mtext("d)", cex=0.8, at=0.8, line=0.2)

par(mar=c(2.1,1.8,1.6,0.1))
bp <- barplot(height=wrack.syncexpl.bandavg, beside=TRUE, ylim=c(-30,110),
              names.arg = c("2-8 mo", "8-16 mo", "16-60 mo"), col=pal)
abline(h=0, col="black")
abline(v=c(6.5, 12.5), col="grey", lty=2)
legend("topleft", fill=pal[1:4], legend=c("Kelp biomass","Waves","Beach width","Interactions")[1:4],
       bty="n")
legend("bottomleft", fill=pal[5], legend=c("Kelp biomass","Waves","Beach width","Interactions","Total")[5],
       bty="n")
text(x=14.5,y=-5,"na",col="grey", pos=3, cex=0.8)
#text(x=2.5,y=-7,"na",col="grey", pos=3, cex=0.8)
mtext("e)", cex=0.8, at=0.35, line=0.2)

mtext("Timescale (months)", 2, outer=TRUE, at=4/5, cex=0.8, line=0.5)
#mtext("Timescale (months)", 1, outer=TRUE, cex=0.8, line=0.5, at=0.5-0.02)
mtext("Annual timescale\nsynchrony",2,outer=TRUE, cex=0.8, line=0, at=5/10)
mtext("Time-averaged\nsynchrony",2,outer=TRUE, at=3/10+0.03, cex=0.8, line=0)
mtext("Synchrony\n contribution (%)",2,outer=TRUE,at=1/10+0.02,cex=0.8,line=0)

image(t(matrix(1:100)), col=colorfill(100), xaxt="n", yaxt="n")
axis(2, at=seq(0,1,length.out=5), labels=round(seq(0,2.18,length.out=5),1), las=2)
image(t(matrix(1:100)), col=colorfill(100), xaxt="n", yaxt="n")
axis(2, at=seq(0,1,length.out=5), labels=round(seq(0,1.6,length.out=5),1), las=2)

dev.off()


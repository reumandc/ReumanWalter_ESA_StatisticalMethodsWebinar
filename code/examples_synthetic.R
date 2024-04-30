#This script implements a bunch of synthetic examples, mostly taken from the
#wsyn vignette, for use in the talk.
#
#Reuman

#***setup

resloc<-"../results/synthetic_results/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***wavelet example

set.seed(101)

#Set up an artificial time series that changes timescale of oscillation in the middle,
#and has some white noise
time1<-1:100; time2<-101:200; times<-c(time1,time2)

ts1<-c(sin(2*pi*time1/15),0*time2) #period 15 initially
ts2<-c(0*time1,sin(2*pi*time2/8)) #then period 8
ts3<-rnorm(200,mean=0,sd=0.5) #add some white noise

ts<-ts1+ts2+ts3 

#Now apply the wavelet transform, obtaining an object of class `wt`. Default parameter
#values for `scale.min`, `scale.max.input`, `sigma` and `f0` are usually good enough for 
#initial data exploration.

ts<-wsyn::cleandat(ts,times,clev=1)$cdat
wtres<-wsyn::wt(ts,times)

#now plot the magnitude and phase and verify that they make sense given how the data were constructed

grDevices::pdf(file=paste0(resloc,"WaveletExample_magnitude.pdf"))
wsyn::plotmag(wtres)
dev.off()

grDevices::pdf(file=paste0(resloc,"WaveletExample_phase.pdf"))
wsyn::plotphase(wtres)
dev.off()



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

#Set up an artificial time series that changes timescale of oscillation in the 
#middle, and has some white noise
time1<-1:100; time2<-101:200; times<-c(time1,time2)

ts1<-c(sin(2*pi*time1/15),0*time2) #period 15 initially
ts2<-c(0*time1,sin(2*pi*time2/8)) #then period 8
ts3<-rnorm(200,mean=0,sd=0.5) #add some white noise

ts<-ts1+ts2+ts3 

#Now apply the wavelet transform, obtaining an object of class `wt`. Default 
#parameter values for `scale.min`, `scale.max.input`, `sigma` and `f0` are usually 
#good enough forinitial data exploration.
ts<-wsyn::cleandat(ts,times,clev=1)$cdat
wtres<-wsyn::wt(ts,times)

#now plot the magnitude and phase and verify that they make sense given how the 
#data were constructed
grDevices::pdf(file=paste0(resloc,"WaveletExample_magnitude.pdf"))
wsyn::plotmag(wtres)
grDevices::dev.off()

grDevices::pdf(file=paste0(resloc,"WaveletExample_phase.pdf"))
wsyn::plotphase(wtres)
grDevices::dev.off()

#***wavelet phasor mean field example

#Set up multi-location data that all has three components: 1) a sin wave of period 
#10 for the first half of the data and then 5 for the second half, present in all\
#locations; 2) a sin wave of period 3, phase randomized independently in all 
#locations; 3) white noise. 
times1<-0:50; times2<-51:100; times<-c(times1,times2)
ts1<-c(sin(2*pi*times1/10),sin(2*pi*times2/5))+1.1 

dat<-matrix(NA,11,length(times))
for (counter in 1:dim(dat)[1])
{
  ts2<-3*sin(2*pi*times/3+2*pi*runif(1))+3.1
  ts3<-rnorm(length(times),0,1.5)
  dat[counter,]<-ts1+ts2+ts3    
}
dat<-wsyn::cleandat(dat,times,1)$cdat

#synchrony is not visually obvious
grDevices::pdf(file=paste0(resloc,"WPMFExample_timeseries.pdf"))
plot(times,dat[1,]/10+1,type='l',xlab="Time",ylab="Time series index",ylim=c(0,12))
for (counter in 2:dim(dat)[1])
{
  lines(times,dat[counter,]/10+counter)
}
grDevices::dev.off()

#Nor can synchrony be readily detected by examining the $55$ pairwise correlation coefficients
cmat<-cor(t(dat))
diag(cmat)<-NA
cmat<-as.vector(cmat)
cmat<-cmat[!is.na(cmat)]
grDevices::pdf(file=paste0(resloc,"WPMFExample_PairwiseCorrelation.pdf"))
hist(cmat,30,xlab="Pearson correlation",ylab="Count")
grDevices::dev.off()

#so now apply the wpmf
wpmfres<-wsyn::wpmf(dat,times,sigmethod="quick")
grDevices::pdf(file=paste0(resloc,"WPMFExample_wpmf.pdf"))
wsyn::plotmag(wpmfres,sigthresh=0.95)
grDevices::dev.off()
#you can see what we know is there because we built it in, so this demonstrates the utility of the method

#***wavelet coherence




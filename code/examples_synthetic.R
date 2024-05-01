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

#make some environmental data, 2 superimposed sin waves in all locations and white noise
times<-(-3:100); lt<-length(times)
ts1<-sin(2*pi*times/10)
ts2<-5*sin(2*pi*times/3)
x<-matrix(ts1+ts2,11,lt,byrow=TRUE)+
  matrix(rnorm(11*lt,0,1.5),11,lt)

#make biological data - the population is a moving average of the environment plus white noise
times<-0:100; lt<-length(times)
y<-matrix(NA,11,lt) #the driven (biological) variable
for (i in 1:101) {
  y[,i]<-apply(FUN=mean,X=x[,i:(i+2)],MARGIN=1) }
y<-y+matrix(rnorm(11*lt,mean=0,sd=3),11,lt)

x<-x[,4:104]
x<-wsyn::cleandat(x,times,1)$cdat 
y<-wsyn::cleandat(y,times,1)$cdat

#The relationship between the environmental and biological variables cannot 
#readily be detected using ordinary correlation methods 
allcors<-c()
for (counter in 1:dim(x)[1])
{
  allcors[counter]<-cor(x[counter,],y[counter,])
}
grDevices::pdf(file=paste0(resloc,"CoherenceExample_LocalBiolEnvCorrs.pdf"))
hist(allcors,xlab="Biol/env correlation",ylab="Count")
grDevices::dev.off()

#However, the function `coh` can be used to compute the coherence, which reveals a 
#timescale-specific relationship
if (file.exists(paste0(resloc,"CoherenceResults.Rds")))
{
  cohres<-readRDS(file=paste0(resloc,"CoherenceResults.Rds"))
} else
{
  cohres<-wsyn::coh(dat1=x,dat2=y,times=times,norm="powall",
         sigmethod="fftsurrog1",nrand=1000,
         f0=0.5,scale.max.input=28)
  saveRDS(cohres,file=paste0(resloc,"CoherenceResults.Rds"))
}

#now display what you get
grDevices::pdf(file=paste0(resloc,"CoherenceExample_PlotCoh1.pdf"))
wsyn::plotmag(cohres)
grDevices::dev.off()

#now add p-values for a couple bands
cohres<-wsyn::bandtest(cohres,c(8,12))
cohres<-wsyn::bandtest(cohres,c(2,4))

#now display what you get again
grDevices::pdf(file=paste0(resloc,"CoherenceExample_PlotCoh2.pdf"))
wsyn::plotmag(cohres)
grDevices::dev.off()

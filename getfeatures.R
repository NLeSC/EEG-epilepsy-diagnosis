getfeatures = function(x,fns) {
  
  #--------------------------------------------------------
  # define generic function for feature extraction
  # Note: approx_entropy from pracma package is very slow
  require(entropy)
  require(signal); require(moments)
  
  # x: dataframe with multiple signals in columns over which the feature is calculated
  # fn: feature names that need to be calculated
  en.entropy = function(x) { # entropy
    en.entropy = entropy::entropy(discretize(as.numeric(x)+runif(n=length(x),min=0,max=1e-8)+0.0001,numBins=100))
  }
  domfreq = function(x) {
    pp = spectrum(x,plot=FALSE)
    domfreq = pp$freq[which.max(pp$spec)]
  }
  maxpow = function(x) {
    pp = spectrum(x,plot=FALSE)
    maxpow = max(pp$spec^2)
  }
  RMS = function(x) {
    RMS = sqrt(mean(x^2))
  }
  zerocross = function(x) { # counts the number of zero crossing
    cc = rep(0,length(x))
    bf.fil = function(y) { # band-pass filter
      sf = 1000
      lb = 0.2
      hb = 60
      Wc = matrix(0,2,1) 
      Wc[1,1] = lb / (sf/2) 
      Wc[2,1] = hb / (sf/2) 
      n = 4
      bf = signal::butter(n,Wc,type=c("pass")) 
      bf.fil = signal::filter(bf,y) 
    }
    x = bf.fil(x)
    cc[which(x >= 0)] = 1
    cc[which(x < 0)] = 0
    zerocross = length(which(abs(diff(cc)) == 1))
  }
  # t = c()
  if (length(which(fns=="mean")) >0) A=aggregate(x=x,by=list(x$bands),mean)
  if (length(which(fns=="sd")) >0) A=aggregate(x=x,by=list(bands),sd)
  if (length(which(fns=="entropy")) >0) {
    A=aggregate(x=x,by=list(bands),en.entropy)
  }
    # if (length(which(fn=="en.entropy")) >0) t = rbind(t,sapply(x,en.entropy))
  #   if (length(which(fn=="max")) >0) t = rbind(t,sapply(x,max))
  #   if (length(which(fn=="min")) >0) t = rbind(t,sapply(x,min))
  #   if (length(which(fn=="skewness")) >0) t = rbind(t,sapply(x,moments::skewness))
  #   if (length(which(fn=="median")) >0) t = rbind(t,sapply(x,median))
  #   if (length(which(fn=="domfreq")) >0) t = rbind(t,sapply(x,domfreq))
  #   if (length(which(fn=="maxpow")) >0) t = rbind(t,sapply(x,maxpow))
  #   if (length(which(fn=="zerocross")) >0) t = rbind(t,sapply(x,zerocross))
  if (length(which(fn=="RMS")) >0) A=aggregate(x=x,by=list(bands),RMS)
  # A = as.data.frame(t)
  return(A)
}
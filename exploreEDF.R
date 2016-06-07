rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
fname = "data/emotivexample/example-data-control--26.05.2016.22.37.50.edf"
require(edfReader)
H = readEdfHeader(fname)
S = readEdfSignals(H,fragments=TRUE)
samplerate = S$AF3$sRate
signal = S$AF3$signal
labels = names(S)
nLabels = length(labels)
cnt = 8
for (i in 3:16) {
  if (cnt == 8) {
    x11()
    par(mfrow=c(4,2))
    cnt = 0
  }
  x = S[[labels[i]]]$signal
  starttime = as.POSIXlt(as.numeric(S[[labels[i]]]$startTime),format="%d-%m-%Y %H:%M:%S",origin="01-01-1966 00:00:00 CEST")
  ts = seq(from=starttime,by=1/128,length.out=length(x))
  plot(ts,x,type="l",main=labels[i])
  cnt = cnt +1
}
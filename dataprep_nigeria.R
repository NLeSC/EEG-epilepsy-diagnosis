rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
source("clean_emotiv.R")
source("extract_features.R")
library(wavelets)
library(pracma)

doclean = TRUE
extractfeature = TRUE
sf = 128 #sample frequency
epochlength = 4 # in seconds
if (doclean == TRUE) {
  datadir =  "/media/windows-share/EEGs_Nigeria" #"data/eeg"
  metadatafile = "/media/windows-share/merged-meta-data_nigeria.csv"
  outputdir =  "/media/windows-share/EEGs_Nigeria_cleaned" 
  gyrothreshold = 30 #gyro unit deviations from the median that are not tolerated
  mindur = 4 # minimum duration in seconds
  # define known errors based on Research Remarks (this is hardcoded, but could be extracted from a file in the future)
  knownerrors = list(c(84,"0:40","0:50"), # => id, starttime, endtime of the problematic period
                     c(91,"0:40","0:50")) # if there are more periods per id then use new entry
  referencegroup = "control"
  condition_start_closed = "closed"
  protocolvariable = "first_condition"
  protocoltimes = c(30,150,270) # in seconds
  amountdata = clean_emotiv(datadir,metadatafile,outputdir,sf,gyrothreshold,mindur,knownerrors,
                            protocoltimes,referencegroup,condition_start_closed)
  print(paste0("successful open: ",length(which(is.na(amountdata[,1]) == FALSE)) / nrow(amountdata)))
  print(paste0("successful closed: ",length(which(is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
  print(paste0("succesful both: ",length(which(is.na(amountdata[,1]) == FALSE &
                                                 is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
}
if (extractfeature == TRUE) {
  datadir = "/media/windows-share/EEGs_Nigeria_cleaned" 
  # featurenames and wavelet filter types to be extracted:
  # fn = c("mean","sd","entropy","max","min","skewness",
  #        "median","domfreq","maxpow","zerocross","RMS","pracma.samen") # feature names ,"lyapunov"
  fn = c("sd","entropy","min","max") #
  # filtertypes =  c(paste0("d",seq(2,20,by=2)), # Daubechies
  #                  paste0("la",seq(8,20,by=2)), #Least Aymetric
  #                  paste0("bl",c(14,18,20)), #Best localized
  #                  paste0("c",seq(6,30,by=6))) # Coiflet
  filtertypes =  c(paste0("d",seq(2,20,by=4)), # Daubechies
                   paste0("la",seq(8,20,by=2)), #Least Aymetric
                   paste0("bl",c(14,18,20)), #Best localized
                   paste0("c",seq(6,30,by=6))) # Coiflet
  # filtertypes =  c("d4") # Daubechies
  n.levels = 7
  ef = extract_features(datadir,sf,n.levels,filtertypes,epochlength)
  DAT = ef$DAT
  LAB = ef$LAB
  save(DAT,LAB,labels,file=paste0("data/features_nigeria_",epochlength,".RData"))
}
rm(list=ls())
graphics.off()
library(wavelets)
library(pracma)
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
funcfiles = list.files("functions",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) {
  source(i)
}


doclean = FALSE
extractfeature = TRUE
sf = 128 #sample frequency
epochlength = 10 # in seconds
if (doclean == TRUE) {
  datadir =  "/media/windows-share/EEGs_Nigeria" #"data/eeg"
  metadatafile = "/media/windows-share/merged-meta-data_nigeria.csv"
  outputdir =  "/media/windows-share/EEGs_Nigeria_cleaned" 
  gyrothreshold = 30 #gyro unit deviations from the median that are not tolerated
  mindur = 4 # minimum duration in minutes
  # define known errors based on Research Remarks (this is hardcoded, but could be extracted from a file in the future)
  knownerrors = list(c(84,"0:30","0:50"), # => id, starttime, endtime of the problematic period
                     c(91,"0:30","0:50"), # if there are more periods per id then use new entry
                     c(104,"0:30","0:55"), c(104,"1:02","1:13"), c(105,"2:30","2:55"), c(108,"2:30","2:52"),
                     c(134,"2:30","2:47"), c(500,"2:30","2:46"), c(508,"0:30","4:30"), c(511,"2:30","2:47"),
                     c(515,"0:30","4:30"),c(519,"2:30","3:00"), c(526,"2:30","4:30"), c(534,"2:30","2:47"),
                     c(569,"2:30","3:10"), c(580,"2:45","4:30"), c(590,"2:30","2:50"), c(599,"2:30","2:50"),
                     c(622,"0:30","2:30"),c(623,"2:30","2:50"), c(626,"2:15","2:30"))
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
  kkk
  DAT = ef$DAT
  LAB = ef$LAB
  save(DAT,LAB,labels,file=paste0("data/features_nigeria_",epochlength,".RData"))
}
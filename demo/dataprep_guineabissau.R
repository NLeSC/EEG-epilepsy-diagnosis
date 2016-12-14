rm(list=ls())
graphics.off()
library(wavelets)
library(pracma)
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
funcfiles = list.files("functions",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) {
  source(i)
}
doclean = FALSE # only needs to be done once, because it is not depenedent on epochlength
extractfeature = TRUE
sf = 128 #sample frequency
epochlength = 10 # in seconds
if (doclean == TRUE) {
  datadir =  "/media/windows-share/EEGs_Guinea-Bissau__16-06-2016" #"data/eeg"
  metadatafile = "/media/windows-share/subject.id_with_meta-info__anonymized.csv"
  outputdir =  "/media/windows-share/EEGs_Guinea-Bissau_cleaned" 
  gyrothreshold = 30 #gyro unit deviations from the median that are not tolerated
  mindur = 4 # minimum duration in minutes
  # define known errors based on Research Remarks (this is hardcoded, but could be extracted from a file in the future)
  knownerrors = list(c(10,"4:30","5:00"), # => id, starttime, endtime of the problematic period
                     c(13,"3:00","3:30"), # if there are more periods per id then use new entry
                     c(20,"0:55","1:10"),
                     c(22,"2:55","3:10"),
                     c(29,"3:00","3:30"),
                     c(40,"3:00","3:30"),
                     c(44,"3:00","4:00"),
                     c(58,"3:00","3:30"),
                     c(61,"3:00","5:00"),
                     c(67,"3:00","3:30"),
                     c(84,"2:20","2:30"),
                     c(92,"2:10","2:20"))
  referencegroup = "Control"
  condition_start_closed = "closed-3min-then-open-2min"
  protocolvariable = "Eyes.condition"
  protocoltimes = c(60,180,300) # in seconds
  amountdata = clean_emotiv(datadir,metadatafile,outputdir,sf,gyrothreshold,mindur,knownerrors,
                            protocoltimes,referencegroup,condition_start_closed,protocolvariable)
  print(paste0("successful open: ",length(which(is.na(amountdata[,1]) == FALSE)) / nrow(amountdata)))
  print(paste0("successful closed: ",length(which(is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
  print(paste0("succesful both: ",length(which(is.na(amountdata[,1]) == FALSE &
                                                 is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
}
if (extractfeature == TRUE) {
  datadir = "/media/windows-share/EEGs_Guinea-Bissau_cleaned"
  # featurenames and wavelet filter types to be extracted:
  # fn = c("mean","sd","entropy","max","min","skewness",
  #        "median","domfreq","maxpow","zerocross","RMS","pracma.samen") # feature names ,"lyapunov"
  fn = c("sd","entropy","min","max") #
  # filtertypes =  c(paste0("d",seq(2,20,by=2)), # Daubechies
  #                  paste0("la",seq(8,20,by=2)), #Least Aymetric
  #                  paste0("bl",c(14,18,20)), #Best localized
  #                  paste0("c",seq(6,30,by=6))) # Coiflet
#   filtertypes =  c(paste0("d",seq(2,20,by=4)), # Daubechies
#                    paste0("la",seq(8,20,by=2)), #Least Aymetric
#                    paste0("bl",c(14,18,20)), #Best localized
#                    paste0("c",seq(6,30,by=6))) # Coiflet
  filtertypes =  paste0("d",seq(2,10,by=4)) # Daubechies
  n.levels = 7
  ef = extract_features(datadir,sf,n.levels,filtertypes,epochlength)
  DAT = ef$DAT
  LAB = ef$LAB
  save(DAT,LAB,labels,file=paste0("features_and_bestmodels/features_ginneabissau_",epochlength,".RData"))
}
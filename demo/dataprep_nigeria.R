rm(list=ls())
graphics.off()
library(wavelets)
library(pracma)
setwd("/home/vincent/utrecht")
shareddrive = "/media/windows-share/EEG"
funcfiles = list.files("emotivepilepsy/R",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) source(i)

doclean = TRUE
extractfeature = FALSE
sf = 128 #sample frequency

if (doclean == TRUE) {
  datadir =  paste0(shareddrive,"/EEGs_Nigeria") #"data/eeg"
  metadatafile = paste0(shareddrive,"/merged-meta-data_nigeria.csv")
  outputdir =  paste0(shareddrive,"/EEGs_Nigeria_cleaned")
  gyrothreshold = 30 #gyro unit deviations from the median that are not tolerated
  mindur = 4 # minimum duration in minutes
  # define known errors based on Research Remarks (this is hardcoded, but could be extracted from a file in the future)
  knownerrors = list(c(10,"2:35","2:45"),
                     # file 11.1 was removed following the remarks by the researchers
                     c(40,"2:35","2:45"),
                     c(63,"2:15","2:25"),
                     # file 65.1, 65.2 , and 69 removed following the remarks by the researchers
                     c(84,"0:30","0:50"), # => id, starttime, endtime of the problematic period
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
  clean_stats = clean_emotiv(datadir,metadatafile,outputdir,sf,gyrothreshold,mindur,knownerrors,
                             protocoltimes,referencegroup,condition_start_closed,protocolvariable)
  amountdata = clean_stats$amountdata
  correction_overview = clean_stats$correction_overview
  save(correction_overview,
       file=paste0(shareddrive,"/features_and_bestmodels/correctionoverview_nigeria.RData"))
  print(paste0("successful open: ",length(which(is.na(amountdata[,1]) == FALSE)) / nrow(amountdata)))
  print(paste0("successful closed: ",length(which(is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
  print(paste0("succesful both: ",length(which(is.na(amountdata[,1]) == FALSE &
                                                 is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
}
for (epochlength in c(10,4)) { # in seconds
  if (extractfeature == TRUE) {
    datadir = paste0(shareddrive,"/EEGs_Nigeria_cleaned")
    # featurenames and wavelet filter types to be extracted:
    fn = c("sd") #,"entropy","min","max") #
    filtertypes =  paste0("d",seq(2,10,by=4)) # Daubechies
    n.levels = 7
    ef = extract_features(datadir,sf,n.levels,filtertypes,epochlength,fn)
    DAT = ef$DAT
    LAB = ef$LAB
    save(DAT,LAB,labels,file=paste0(shareddrive,"/features_and_bestmodels/features_nigeria_",epochlength,".RData"))
  }
}
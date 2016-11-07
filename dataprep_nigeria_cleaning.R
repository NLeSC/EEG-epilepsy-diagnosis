rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
datadir =  "/media/windows-share/EEGs_Nigeria" #"data/eeg"
metadatafile = "/media/windows-share/merged-meta-data_nigeria.csv"
outputdir =  "/media/windows-share/EEGs_Nigeria_cleaned" 
sf = 128 #sample frequency
gyrothreshold = 30 #gyro unit deviations from the median that are not tolerated
mindur = 4 # minimum duration in seconds

source("clean_emotiv.R")

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
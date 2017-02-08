#==============================================
# Update the following lines:

setwd("/home/vincent/utrecht")
knownerrorfile = "/media/windows-share/EEG/nigeria_knownerrors.csv"
metadatafile = "/media/windows-share/EEG/merged-meta-data_nigeria.csv"
datadir = "/media/windows-share/EEG/EEGs_Nigeria"
outputdir = "/media/windows-share/EEG" # this is where folders will be created to store the output
namecountry = "ni"

referencegroup = "control"
condition_start_closed = "closed"
protocolvariable = "first_condition"
protocoltimes = c(30,150,270) # in seconds
#===============================================
funcfiles = list.files("EEG-epilepsy-diagnosis/R",include.dirs=TRUE,full.names = TRUE) # this line only needed when developing
for (i in funcfiles) source(i) # this line only needed when developing

doclean = TRUE
extractfeature = FALSE

outputdir_clean = paste0(outputdir,"/EEGs_",namecountry,"_cleaned")
if (!file.exists(outputdir_clean))  dir.create(outputdir_clean)
outputdir_features = paste0(outputdir,"/EEGs_features")
if (!file.exists(outputdir_features)) dir.create(outputdir_features)
outputdir_logs = paste0(outputdir,"/EEGs_logs")
if (!file.exists(outputdir_logs)) dir.create(outputdir_logs)
outputdir_images = paste0(outputdir,"/images")
if (!file.exists(outputdir_images)) dir.create(outputdir_images)

if (doclean == TRUE) {
  knownerrors = read.csv(knownerrorfile)  
  clean_stats = clean_emotiv(datadir,metadatafile,outputdir,knownerrors,
                             protocoltimes,referencegroup,condition_start_closed,protocolvariable)
  amountdata = clean_stats$amountdata
  correction_overview = clean_stats$correction_overview
  save(correction_overview,
       file=paste0(outputdir_logs,"/correctionoverview_",nmaecountry,".RData"))
  print(paste0("successful open: ",length(which(is.na(amountdata[,1]) == FALSE)) / nrow(amountdata)))
  print(paste0("successful closed: ",length(which(is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
  print(paste0("succesful both: ",length(which(is.na(amountdata[,1]) == FALSE &
                                                 is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
}
if (extractfeature == TRUE) {
  for (epochlength in c(4)) { # in seconds
    # featurenames and wavelet filter types to be extracted:
    fn = c("sd") #,"entropy","min","max") #
    filtertypes =  paste0("d",seq(2,10,by=4)) # Daubechies
    ef = extract_features(cleandatadir=outputdir_clean,filtertypes,epochlength,fn)
    DAT = ef$DAT
    LAB = ef$LAB
    save(DAT,LAB,labels,file=paste0(outputdir_features,"/features_",namecountr,"_epoch",epochlength,".RData"))
  }
}

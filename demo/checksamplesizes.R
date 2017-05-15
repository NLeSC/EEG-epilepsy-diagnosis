rm(list=ls())
shareddrive = "/media/windows-share/EEG"
nimeta = read.csv(paste0(shareddrive,"/merged-meta-data_nigeria.csv"))
gbmeta = read.csv(paste0(shareddrive,"/subject.id_with_meta-info__anonymized.csv"))

# Error message in GB data:
Error in ROCR::prediction(y_pred, y_true) : 
  Number of classes is not equal to 2.
ROCR currently supports only evaluation of binary classification tasks.

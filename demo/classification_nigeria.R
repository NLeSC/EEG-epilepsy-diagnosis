rm(list=ls())
graphics.off()
# library(caret)
# library(psych)
# library(pROC)
library(emotivepilepsy)
setwd("/home/vincent/utrecht")
# load(file="data/features_nigeria_4.RData"); logdur = 4
load(file="features_and_bestmodels/features_nigeria_10.RData"); logdur = 10
funcfiles = list.files("emotivepilepsy/R",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) source(i)


LAB$diagnosis = as.character(LAB$diagnosis)
LAB$diagnosis[which(LAB$diagnosis == "control")] = "Control"
LAB$diagnosis[which(LAB$diagnosis == "epilepsy")] = "Epilepsy"
LAB$diagnosis = as.factor(LAB$diagnosis)

# proto_i = "eyesclosed"
proto_i = 1 #"eyesopen" # "open" 1 #closed 2
logfile = "data/log_guinneabissau.csv" # not used when uselog = FALSE

perid = FALSE
trainbestmodel = FALSE #option to turn this off for Nigeria

# tidy up formatting to be suitable for classifier training
RDL = reformat_DATLAB(DAT,LAB,aggregateperid=perid) # aggregate per unique id
DAT =RDL$DAT
LAB = RDL$LAB

#===============================================================
# split data in training, validation and test set
P = split_data(LAB,DAT,logfile = logfile,proto_i=proto_i,split=c(20,20),
               uselog = FALSE,logdur=logdur)
LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
# generate dictionary of model characteristics
modeldict = create_modeldict(DAT)
#===============================================================
# train models or loaded previously trained model
if (trainbestmodel == TRUE) {
  trainingresults = train_model(DATtrain,LABtrain,DATval,LABval,modeldict,classifier = "rf",perid=perid)
  best_model = trainingresults$best_model
  modelcomparison = trainingresults$result
  fes = trainingresults$fes
  country = "ni"
  bestmodelfile = paste0("features_and_bestmodels/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",perid,".RData")
  save(best_model,fes,file=bestmodelfile) # Save best model
  rm(best_model,fes)
} else {
  country = "gb" # load model from other country
  bestmodelfile = paste0("features_and_bestmodels/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",perid,".RData")
  LABtest = rbind(LABval) # ignore test data, and only evaluate on training and validation data
  DATtest = rbind(DATval)
}
# Reload best model
load(bestmodelfile)
#===============================================================
# evaluate on test set
test_factors = DATtest[,fes]
evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=perid)
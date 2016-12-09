rm(list=ls())
graphics.off()
library(caret)
library(psych)
library(pROC)

setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
# load(file="data/features_ginneabissau_4.RData"); logdur = 4
load(file="data/features_ginneabissau_10.RData"); logdur = 10
funcfiles = list.files("functions",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) {
  source(i)
}
# proto_i = "eyesclosed"
proto_i = 2 #"eyesopen" #"open" =1 #closed= 2
logfile = "data/log_guinneabissau.csv" # not used when uselog = FALSE

perid = TRUE
trainbestmodel = TRUE #option to turn this off for Nigeria

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
# train and evaluate all models
if (trainbestmodel == TRUE) {
  trainingresults = train_model(DATtrain,LABtrain,DATval,LABval,modeldict)
  best_model = trainingresults$best_model_randomforest
  modelcomparison = trainingresults$result
  fes = trainingresults$fes
  country = "gb" #"gb" =  guinea bisau
  bestmodelfile = paste0("data/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",perid,".RData")
  save(best_model,fes,file=bestmodelfile) # Save best model
  rm(best_model,fes)
} else {
  country = "ni" # load model from other country
  bestmodelfile = paste0("data/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",perid,".RData")
    LABtest = rbind(LABtrain,LABval) # ignore test data, and only evaluate on training and validation data
    DATtest = rbind(DATtrain,DATval)
#   reshuffle = sample(1:nrow(DATval))
#   DATtest = DATtest[reshuffle,]
#   LABtest = LABtest[reshuffle,]
}
# Reload best model
load(bestmodelfile)
#===============================================================
# evaluate on test set
test_factors = DATtest[,fes]
pred_test = predict(best_model, test_factors,type="prob")
result.roc <- roc(DATtest$diagn, pred_test$X1)
# Code to plot the ROC curve if we wanted to:
# x11()
# plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auctest = result.roc$auc
result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
pred_test_cat = rep("X1",nrow(pred_test))
pred_test_cat[pred_test$X2 > 0.500] = "X2"
confmat = create_confmatrix(pred_test_cat,LABtest$diagnosis) 
test.confmatrix = paste0(confmat[1,1:2],"_",confmat[2,1:2],collapse="_")
test.auc = round(auctest,digits=3)
test.kappa = round(cohen.kappa(x=confmat)$kappa,digits=3)
test.acc = round(sum(diag(confmat)) / sum(confmat),digits=3)
print(paste0(proto_i," acc ",test.acc," kappa ",test.kappa," auc ",test.auc," ",test.confmatrix))
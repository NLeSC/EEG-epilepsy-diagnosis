rm(list=ls())
graphics.off()
library(caret)
library(psych)
# library(ROCR)
library(pROC)

setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
# load(file="data/features_ginneabissau_4.RData"); logdur = 4
load(file="data/features_ginneabissau_10.RData"); logdur = 10
funcfiles = list.files("functions",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) {
  source(i)
}

# proto_i = "eyesclosed"
proto_i = "eyesopen" #"open" #closed
logfile = "data/log_guinneabissau.csv" # not used when uselog = FALSE

#===============================================================
# split data in training, validation and test set
P = split_data(LAB,DAT,logfile = logfile,proto_i=proto_i,split=c(20,20),uselog = FALSE,logdur=logdur,useallepoch=TRUE)
LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
kkk
# generate dictionary of model characteristics
modeldict = create_modeldict(DAT)
#===============================================================
# train and evaluate all models
trainingresults = train_model(DATtrain,LABtrain,DATval,LABval,modeldict)
best_model = trainingresults$best_model_randomforest
modelcomparison = trainingresults$result
fes = trainingresults$fes


country = "gb"
bestmodelfile = paste0("data/bestmodel_",proto_i,"_dur",logdur,"_country",country,".RData")
# Save best model
save(best_model,fes,file=bestmodelfile)
rm(best_model,fes)
# Reload best model
load(bestmodelfile)

#===============================================================
# evaluate on test set
test_factors = DATtest[,fes]
pred_test = predict(best_model, test_factors,type="prob")
result.roc <- roc(DATtest$diagn, pred_test$Control)
# Code to plot the ROC curve if we wanted to:
# x11()
# plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auctest = result.roc$auc
result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
pred_test_cat = rep("Control",nrow(pred_test))
pred_test_cat[pred_test$Epilepsy > 0.500] = "Epilepsy"
confmat = create_confmatrix(pred_test_cat,LABtest$diagnosis) 
test.confmatrix = paste0(confmat[1,1:2],"_",confmat[2,1:2],collapse="_")
test.auc = round(auctest,digits=3)
test.kappa = round(cohen.kappa(x=confmat)$kappa,digits=3)
test.acc = round(sum(diag(confmat)) / sum(confmat),digits=3)
print(paste0(proto_i," acc ",test.acc," kappa ",test.kappa," auc ",test.auc," ",test.confmatrix))
#======================================
# Results:


# eyes closed, 10 second epoch
# 1epoch per person: "eyesclosed acc 0.7 kappa 0.4 auc 0.77 7_3_3_7"

# eyes open, 10 second epoch
# "eyesopen acc 0.8 kappa 0.6 auc 0.91 9_1_3_7"



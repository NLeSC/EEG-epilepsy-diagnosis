rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
load(file="data/features_ginneabissau.RData")

LOG = read.csv("log.csv",stringsAsFactors = FALSE)

# make sure that both LAB and DAT have matching row order
getid = function(x) {
  tmp = unlist(strsplit(x,"_"))[2]
  tmp2 = unlist(strsplit(tmp,"id"))[2]
  return(as.numeric(tmp2))
}
DAT$id = as.numeric(sapply(DAT$fnames,getid))
DAT = DAT[order(DAT$id),]
if (length(which((rownames(LAB) == DAT$fnames) == FALSE)) > 0) print("Error: data order does not match")

# DAT = DATmax[,c(seq(1,50,by=5),seq(4,50,by=5),seq(5,50,by=5))]
# merge min and max features
DAT = DAT[,-c(ncol(DAT),ncol(DAT)-1)]

# investigate correlation between two implementations of entropy 
# CMAX = rep(0,7)
# CMAX[1] = cor(DAT$'1.entropy.d4.max',DAT$'1.pracma.samen.d4.max')
# CMAX[2] = cor(DAT$'2.entropy.d4.max',DAT$'2.pracma.samen.d4.max')
# CMAX[3] = cor(DAT$'3.entropy.d4.max',DAT$'3.pracma.samen.d4.max')
# CMAX[4] = cor(DAT$'4.entropy.d4.max',DAT$'4.pracma.samen.d4.max')
# CMAX[5] = cor(DAT$'5.entropy.d4.max',DAT$'5.pracma.samen.d4.max')
# CMAX[6] = cor(DAT$'6.entropy.d4.max',DAT$'6.pracma.samen.d4.max')
# CMAX[7] = cor(DAT$'7.entropy.d4.max',DAT$'7.pracma.samen.d4.max')

varnames = names(DAT)

getfeaturetype = function(x) {
  return(unlist(strsplit(x,"[.]"))[2])
}
getwavelettype = function(x) {
  tt = unlist(strsplit(x,"[.]"))
  return(tt[length(tt)-1])
}
v2 = sapply(varnames,FUN = getwavelettype)
v3 = sapply(varnames,FUN = getfeaturetype)

lookatprotocol = "closed"
# training set
train_fnames = LOG[which(LOG$set == "train" & LOG$dur == 10 & LOG$protocol == lookatprotocol),]$filename
traini = which((rownames(LAB) %in% train_fnames) == TRUE)
LABtrain = LAB[traini,]
DATtrain = DAT[traini,]
# test set
test_fnames = LOG[which(LOG$set == "test" & LOG$dur == 10 & LOG$protocol == lookatprotocol),]$filename
testi = which((rownames(LAB) %in% test_fnames) == TRUE)
LABtest = LAB[testi,]
DATtest = DAT[testi,]
# validation set
val_fnames = LOG[which(LOG$set == "valid" & LOG$dur == 10 & LOG$protocol == lookatprotocol),]$filename
vali = which((rownames(LAB) %in% val_fnames) == TRUE)
LABval = LAB[vali,]
DATval = DAT[vali,]
# define models to be trained
uv2 = unique(v2)
uv3 = unique(v3)
uv4 = c("rf","glm")
Ncases = length(uv2) #* length(uv4) #* length(uv3) 
result = data.frame(wavelet=rep(" ",Ncases),
                    feature=rep(" ",Ncases),
                    model=rep(" ",Ncases),
                    stringsAsFactors=FALSE)
# DAToriginal = DAT
DATtest$diagn = as.factor(LABtest$diagnosis)
DATtrain$diagn = as.factor(LABtrain$diagnosis)
DATval$diagn = as.factor(LABval$diagnosis)
cnt = 1
library(caret)
require(psych)
# library(ROCR)
library(pROC)
for (wtype in uv2) {
  # for (ftype in uv3) {
  
  ctrl = trainControl(method = "none") #"repeatedcv",number=3,repeats=1
  # ,  savePredictions = TRUE,classProbs=TRUE,summaryFunction = twoClassSummary
  for (modeli in c("rf")) { #,"glm"
    print(paste0(wtype," ",modeli))
    #============================================================
    # training
    if (modeli == "rf") {
      # random forest
      grid_rf = expand.grid(.mtry = c(2)) #2,4,8,16
      m_rf = train(diagn ~ .,data=DATtrain[,c(which(v2==wtype),ncol(DATtrain))],method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)
    } else {
      m_rf = train(diagn ~ .,  data=DATtrain[,c(which(v2==wtype),ncol(DATtrain))], method="glm", metric="Kappa",family="binomial")
    }
    # result$training.kappa[cnt] = round(max(m_rf$results$Kappa),digits=3)
    #===========================================================
    # TO DO: implement hyperparameter optimization using the validation set
    # testing
    pred_val = predict(m_rf, DATval[,1:(ncol(DATtest)-1)],type="prob")
    result.roc <- roc(DATval$diagn, pred_val$Control)
    aucval = result.roc$auc
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
    result$val.auc[cnt] = round(aucval,digits=3)
    pred_val_cat = rep("Control",nrow(pred_val))
    pred_val_cat[pred_val$Epilepsy > 0.500] = "Epilepsy"
    confmat = table(pred_val_cat,LABval$diagnosis)
    result$val.kappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=3)
    #============================================================
    # testing
    pred_test = predict(m_rf, DATtest[,1:(ncol(DATtest)-1)],type="prob")
    result.roc <- roc(DATtest$diagn, pred_test$Control)
    # Code to plot the ROC curve if we wanted to:
    #         x11()
    #         plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
    auctest = result.roc$auc
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
    result$test.auc[cnt] = round(auctest,digits=3)
    # Kappa statistic as an alternative performance metric:
    pred_test_cat = rep("Control",nrow(pred_test))
    pred_test_cat[pred_test$Epilepsy > 0.500] = "Epilepsy"
    confmat = table(pred_test_cat,LABtest$diagnosis)
    result$test.kappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=3)
    #===================================================================
    # keep track of wavelettype and model
    result$wavelet[cnt] = wtype
    result$model[cnt] = modeli
    print(paste0(modeli," ",wtype," validation K:",result$val.kappa[cnt],
                 " test K:",result$test.kappa[cnt],
                 " test auc:",result$test.auc[cnt]))
    # jjj
    cnt = cnt+1
  }
}
result = result[with(result,order(val.kappa)),]
# result$kappa = round(result$kappa,digits=3)
# result$trainingkappa  = round(result$trainingkappa,digits=3)
write.csv(result,file="data/result_guineabissau.csv")

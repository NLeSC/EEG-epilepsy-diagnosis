rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
load(file="data/features.RData")
# DAT = DATmax[,c(seq(1,50,by=5),seq(4,50,by=5),seq(5,50,by=5))]
# merge min and max features
varnames = names(DAT)

getfeaturetype = function(x) {
  return(unlist(strsplit(x,"[.]"))[2])
}
getwavelettype = function(x) {
  return(unlist(strsplit(x,"[.]"))[3])
}
v2 = sapply(varnames,FUN = getwavelettype)
v3 = sapply(varnames,FUN = getfeaturetype)

cut = which(LAB$definitive_diagnosis == 3)
DAT = DAT[-cut,]
LAB = LAB[-cut,]

#define subset:
noepi = which(LAB$definitive_diagnosis == 1)
epi = which(LAB$definitive_diagnosis == 2)
all = c(noepi,epi)
set.seed(300)
traini = c(sample(x=noepi,size=round(length(noepi)*0.5)),
           sample(x=epi,size=round(length(epi)*0.5)))
A =  1:length(LAB$definitive_diagnosis)
testi = which(A %in% traini == FALSE)
LAB$definitive_diagnosis[which(LAB$definitive_diagnosis == 1)] = "NO"
LAB$definitive_diagnosis[which(LAB$definitive_diagnosis == 2)] = "YES"
LABtest = LAB[testi,]
LAB = LAB[traini,]
uv2 = unique(v2)
uv3 = unique(v3)
uv4 = c("rf","glm")
Ncases = length(uv2) * length(uv4) #* length(uv3) 
result = data.frame(wavelet=rep(" ",Ncases),
                    feature=rep(" ",Ncases),
                    model=rep(" ",Ncases),
                    stringsAsFactors=FALSE)
DAToriginal = DAT

cnt = 1
library(caret)
require(psych)
for (wtype in uv2) {
  # for (ftype in uv3) {
    DAT = DAToriginal[,which(v2 == wtype )] #& v3 == ftype
    # take subset and add labels
    DATtest = DAT[testi,]
    DAT = DAT[traini,]
    
    #==========================================
    # Consider starting with PCA to reduce size of data:
#     PCA = prcomp(x=t(DAT),center=TRUE,scale=TRUE,na.action=na.omit, retx=TRUE,tol=0.05)
#     tmp = summary(PCA)
#     PCArot = as.data.frame(PCA$rotation)
#     nPC = which(summary(PCA)$importance[3,] > 0.90)[1]
#     print(nPC)
#     DAT = PCArot[,1:nPC]
    
    # library(pROC)
    #========================================
    
    
    DATtest$diagn = as.factor(LABtest$definitive_diagnosis)
    DAT$diagn = as.factor(LAB$definitive_diagnosis)

    # Classification, caret
#     ctrl = trainControl(method = "repeatedcv",number=5,repeats=2,savePrediction = T,classProbs=TRUE,
#                         summaryFunction=twoClassSummary) #
    ctrl = trainControl(method = "repeatedcv",number=5,repeats=1) #
                        
    for (modeli in c("rf","glm")) {
      if (modeli == "rf") {
        # random forest
        grid_rf = expand.grid(.mtry = c(2,4,8,16))
        m_rf = train(diagn ~ .,data=DAT,method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)
      } else {
        m_rf = train(diagn ~ .,  data=DAT, method="glm", metric="Kappa",family="binomial")
      }
      m_rf_test = predict(m_rf, DATtest)
      confmat = table(m_rf_test,LABtest$definitive_diagnosis)
      #============================================================
      result$trainingkappa[cnt] = round(max(m_rf$results$Kappa),digits=2)
      result$testkappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=2)
      result$wavelet[cnt] = wtype
      # result$feature[cnt] =  ftype
      result$model[cnt] = modeli
      # print(paste0(modeli," ",ftype," ",wtype," training:",result$trainingkappa[cnt]," test:",result$testkappa[cnt]))
      print(paste0(modeli," ",wtype," training:",result$trainingkappa[cnt]," test:",result$testkappa[cnt]))
      cnt = cnt+1
      # print(twoClassSummary(m_rf$pred,lev=levels(m_rf$pred$obs))) # now get kappa
      
#       selectedI = which(m_rf$pred$mtry == 16)
#       plot.roc(as.numeric(m_rf$pred$obs[selectedI]),as.numeric(m_rf$pred$pred[selectedI]))
# #       print(twoClassSummary(m_rf$pred,lev=levels(m_rf$pred$obs)))
#       kjj
    }
  }
# }
result = result[with(result,order(trainingkappa)),]
# result$kappa = round(result$kappa,digits=3)
# result$trainingkappa  = round(result$trainingkappa,digits=3)

write.csv(result,file="data/result.csv")

# # nnet
# grid_nnet = expand.grid(.size = c(4,6),.decay=c(5e-4,5e-5)) 
# m_nnet = train(diagn ~ .,data=DAT,method="nnet",metric="Kappa",trControl=ctrl,tuneGrid=grid_nnet)
# print(m_nnet)
# # C5.0
# grid_c50 = expand.grid(.model = "tree",.trials = c(10,20,30,40),.winnow= "FALSE")
# m_c50 = train(diagn ~ .,data=DAT,method="C5.0",metric="Kappa",trControl=ctrl,tuneGrid=grid_c50)
# print(m_c50)
# #mlp: Kappa = 0 ?
# grid_mlp = expand.grid(.size=4)
# m_mlp = train(diagn ~ .,data=DAT,method="mlp",metric="Kappa",trControl=ctrl,tuneGrid=grid_mlp)
# print(m_mlp)

# # Neural Network: Kappa = 0?
# grid_mlpML = expand.grid(.layer1=c(2,3,4),.layer2=c(2,3),.layer3=c(2,3))
# m_mlpML = train(diagn ~ .,data=DAT,method="mlpML",metric="Kappa",trControl=ctrl,tuneGrid=grid_mlpML)
# print(m_mlpML)

#svm: Kappa = NA?
# grid_lssvmPoly = expand.grid(.degree=c(2,3),.scale=2)
# m_lssvmPoly = caret::train(diagn ~ .,data=DAT,method="lssvmPoly",metric="Kappa",trControl=ctrl,tuneGrid=grid_lssvmPoly)
# print(m_lssvmPoly)
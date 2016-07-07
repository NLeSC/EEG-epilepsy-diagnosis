rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
load(file="data/features_ginneabissau.RData")
# DAT = DATmax[,c(seq(1,50,by=5),seq(4,50,by=5),seq(5,50,by=5))]
# merge min and max features
DAT = DAT[,-ncol(DAT)]
varnames = names(DAT)

getfeaturetype = function(x) {
  return(unlist(strsplit(x,"[.]"))[2])
}
getwavelettype = function(x) {
  return(unlist(strsplit(x,"[.]"))[3])
}
v2 = sapply(varnames,FUN = getwavelettype)
v3 = sapply(varnames,FUN = getfeaturetype)

cut = which(LAB$protocol == "eyesopen")
DAT = DAT[-cut,]
LAB = LAB[-cut,]


#define subset:
noepi = which(as.character(LAB$diagnosis) == "Control")
epi = which(as.character(LAB$diagnosis) == "Epilepsy")
all = c(noepi,epi)
set.seed(300)
traini = c(sample(x=noepi,size=round(length(noepi)*0.5)),
           sample(x=epi,size=round(length(epi)*0.5)))
A =  1:length(LAB$diagnosis)
testi = which(A %in% traini == FALSE)
# LAB$diagnosis[which(as.character(LAB$diagnosis) == "Control")] = "NO"
# LAB$diagnosis[which(as.character(LAB$diagnosis) == "Epilepsy")] = "YES"
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
  DATtest$diagn = as.factor(LABtest$diagnosis)
  DAT$diagn = as.factor(LAB$diagnosis)
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
    confmat = table(m_rf_test,LABtest$diagnosis)
    #============================================================
    result$trainingkappa[cnt] = round(max(m_rf$results$Kappa),digits=2)
    result$testkappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=2)
    result$wavelet[cnt] = wtype
    result$model[cnt] = modeli
    print(paste0(modeli," ",wtype," training:",result$trainingkappa[cnt]," test:",result$testkappa[cnt]))
    cnt = cnt+1
  }
}
result = result[with(result,order(trainingkappa)),]
# result$kappa = round(result$kappa,digits=3)
# result$trainingkappa  = round(result$trainingkappa,digits=3)
write.csv(result,file="data/result_guineabissau.csv")

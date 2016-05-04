rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
load(file="data/features.RData")
# jjjj
# DAT = DATmax[,c(seq(1,50,by=5),seq(4,50,by=5),seq(5,50,by=5))]
# merge min and max features


# TO DO:
# - implement systematic comparison of different wavelets

#==========================================
# Consider starting with PCA to reduce size of data:
# PCA = prcomp(x=t(DAT),center=TRUE,scale=TRUE,na.action=na.omit, retx=TRUE,tol=0.05)
# tmp = summary(PCA)
# PCArot = as.data.frame(PCA$rotation)
# nPC = which(summary(PCA)$importance[3,] > 0.99)[1]
# print(nPC)
# DAT = PCArot[,1:nPC]

#========================================
# Classification, randomForest
noepi = which(LAB$definitive_diagnosis == 1)
epi = which(LAB$definitive_diagnosis == 2)
all = c(noepi,epi)
set.seed(300)
traini = c(sample(x=noepi,size=round(length(noepi)*0.5)),
           sample(x=epi,size=round(length(epi)*0.5)))
# testi = which(traini %in% 1:length(LAB) == FALSE)

# colnames(DATsd) = paste0(colnames(DATsd),"_sd")
# colnames(DATmed) = paste0(colnames(DATmed),"_med")
# DAT = cbind(DATsd[traini,],DATmed[traini,])
# jjjj
DAT = DAT[traini,]
LAB = LAB[traini,]

i1 = which(LAB$definitive_diagnosis == 1)
i2 = which(LAB$definitive_diagnosis == 2)
LAB$definitive_diagnosis[i1] = paste0("c",LAB$definitive_diagnosis[i1])
LAB$definitive_diagnosis[i2] = paste0("c",LAB$definitive_diagnosis[i2])

DAT$diagn = as.factor(LAB$definitive_diagnosis)

#========================================
# Classification, caret
library(caret)
ctrl = trainControl(method = "repeatedcv",number=10,repeats=1)

# random forest
grid_rf = expand.grid(.mtry = c(2,4,8,16))
m_rf = train(diagn ~ .,data=DAT,method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)
print(m_rf)
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
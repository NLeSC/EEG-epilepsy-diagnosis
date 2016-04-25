rm(list=ls())
graphics.off()
load(file="features.RData")

#==========================================
# Consider starting with PCA to reduce size of data:
PCA = prcomp(x=t(DAT),center=TRUE,scale=TRUE,na.action=na.omit, retx=TRUE,tol=0.05)
tmp = summary(PCA)
PCArot = as.data.frame(PCA$rotation)
nPC = which(summary(PCA)$importance[3,] > 0.9)[1]
print(nPC)
DAT = PCArot[,1:nPC]
#========================================
# Classification, randomForest
noepi = which(LAB$definitive_diagnosis == 1)
epi = which(LAB$definitive_diagnosis == 2)
all = c(noepi,epi)
set.seed(300)
traini = c(sample(x=noepi,size=round(length(noepi)*0.5)),
           sample(x=epi,size=round(length(epi)*0.5)))

# colnames(DATsd) = paste0(colnames(DATsd),"_sd")
# colnames(DATmed) = paste0(colnames(DATmed),"_med")
# DAT = cbind(DATsd[traini,],DATmed[traini,])

DAT = DAT[traini,]
LAB = LAB[traini,]
#========================================
# Classification, caret
library(caret)
ctrl = trainControl(method = "repeatedcv",number=5,repeats=5)
grid_rf = expand.grid(.mtry = c(4,8,16,32,64))
DAT$diagn = as.factor(LAB$definitive_diagnosis)
m_rf = train(diagn ~ .,data=DAT,method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)
print(m_rf)
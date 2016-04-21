rm(list=ls())
graphics.off()
load(file="features.RData")

#==========================================
# Consider starting with PCA to reduce size of data
#...
# PCA = prcomp(x=t(DAT),center=TRUE,scale=TRUE,na.action=na.omit, retx=TRUE,tol=0.01)
# tmp = summary(PCA)
# DAT = as.data.frame(PCA$rotation)
#========================================
# Classification, randomForest
noepi = which(LAB$definitive_diagnosis == 1)
epi = which(LAB$definitive_diagnosis == 2)
all = c(noepi,epi)

traini = c(sample(x=noepi,size=round(length(noepi)*0.4)),
           sample(x=epi,size=round(length(epi)*0.4)))
testi = all[which(all %in% traini == FALSE)]

DAT = DAT[traini,]
LAB = LAB[traini,]
# 
# m  = randomForest(DAT[traini,],LAB$definitive_diagnosis[traini],ntree=500)
# p = predict(m,DAT[testi,],type="response")
# x11();plot(p,LAB$definitive_diagnosis[testi],type="p")
# # # confusion matrix
# print(table(round(p),LAB$definitive_diagnosis[testi]))
#========================================
# Classification, caret
print("classification")
library(caret)
ctrl = trainControl(method = "repeatedcv",number=10,repeats=10)
grid_rf = expand.grid(.mtry = c(4,8,16,32))
# DAT = cbind(DAT,)

if (length(which(LAB$definitive_diagnosis == 3)) > 0) { #ignore class 3
  DAT = DAT[-which(LAB$definitive_diagnosis == 3),]
  LAB = LAB[-which(LAB$definitive_diagnosis == 3),]
}
DAT$diagn = as.factor(LAB$definitive_diagnosis)
set.seed(300)
m_rf = train(diagn ~ .,data=DAT,method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)

print(m_rf)
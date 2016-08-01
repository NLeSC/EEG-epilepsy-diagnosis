rm(list=ls())
graphics.off()
library(caret)
library(psych)
# library(ROCR)
library(pROC)

setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
# load(file="data/features_ginneabissau_4.RData")
load(file="data/features_ginneabissau_10.RData")

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

makeconfmatsquare = function(x) { # function to turn confusion matrix square if needed
  if (ncol(x) > nrow(x)) {
    ihave = which(colnames(x) == rownames(x))
    if (ihave == 2) {
      x = rbind(c(0,0),x)
    } else {
      x = rbind(x,c(0,0))
    }
  } else if (ncol(x) < nrow(x)) {
    ihave = which(rownames(x) == colnames(x))
    if (ihave == 2) {
      x = cbind(c(0,0),x)
    } else {
      x = cbind(x,c(0,0))
    }
  }
  return(x)
}
getfeaturetype = function(x) return(unlist(strsplit(x,"[.]"))[2])
getwavelettype = function(x) {
  tt = unlist(strsplit(x,"[.]"))
  return(tt[length(tt)-1])
}
getwaveletlevel= function(x) {
  tt = unlist(strsplit(x,"[.]"))
  return(tt[1])
}
getaggtype= function(x) {
  tt = unlist(strsplit(x,"[.]"))
  return(tt[4])
}

modeldict = data.frame(v2 = sapply(varnames,FUN = getwavelettype),
                       v3 = sapply(varnames,FUN = getfeaturetype),
                       v5 = sapply(varnames,FUN = getwaveletlevel),
                       v6 = sapply(varnames,FUN = getaggtype),stringsAsFactors=FALSE)
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
DATtest$diagn = as.factor(LABtest$diagnosis)
DATtrain$diagn = as.factor(LABtrain$diagnosis)
DATval$diagn = as.factor(LABval$diagnosis)

usePCA = FALSE

# First select best possible features set
testpart = c("wavelet") #,"features","aggregationtype") #,"waveletlevel"
for (testparti in testpart) {
  cat("----------------\n")
  cat(paste0(testparti,"\n"))
  cnt = 1
  #============================================================
  # training
  set.seed(300)
  if (testparti == "wavelet") {
    allvalues = modeldict$v2
  } else if (testparti == "features") {
    allvalues = modeldict$v3
  } else if (testparti == "waveletlevel") {
    allvalues = modeldict$v5
  } else if (testparti == "aggregationtype") {
    allvalues = modeldict$v6
  }
  valueseval = unique(allvalues)
  Ncases = length(valueseval)
  result = data.frame(model=rep(" ",Ncases),
                      stringsAsFactors=FALSE)
  for (jj in 1:length(valueseval)) {
    result$model[cnt] = valueseval[jj]
    fes = which(allvalues==valueseval[jj])
    if (usePCA == TRUE) { # maybe need to replace this with lasso PCA??
      #run pca on training data
      pcatr = prcomp(x=DATtrain[,fes],center=TRUE,scale=TRUE, retx=TRUE,tol=0.05) #,na.action=na.omit
      # now trim it to only have PCAs that explain > 80% of the data
      cumimportance = summary(pcatr)$importance[3,]
      cut = which(cumimportance > 0.7)[1]
      if (cut == 1) cut = 2
      cut = cut - 1
      Npca = cut
      # get pca for each set:
      train_factors = matrix(0,nrow(pcatr$x),Npca)
      val_factors = matrix(0,nrow(DATval[,fes]),Npca)
      test_factors = matrix(0,nrow(DATval[,fes]),Npca)
      for (gi in 1:Npca) {
        train_factors[,gi] = scale(as.matrix(DATtrain[,fes]),center = pcatr$center,scale=pcatr$scale) %*% pcatr$rotation[,gi]
        val_factors[,gi] = scale(as.matrix(DATval[,fes]),center = pcatr$center,scale=pcatr$scale) %*% pcatr$rotation[,gi]
        test_factors[,gi] = scale(as.matrix(DATtest[,fes]),center = pcatr$center,scale=pcatr$scale) %*% pcatr$rotation[,gi]
      }
    } else {
      train_factors = DATtrain[,fes]
      val_factors = DATval[,fes]
      test_factors = DATtest[,fes]
    }
    #========================================================
    # now use factors for random forest
    # set training paramerters
    ctrl = trainControl(method = "repeatedcv",number=10,repeats=3,search="random")
    # train on training set
    m_rf = train(y=DATtrain$diagn,x=train_factors,
                 method="rf",metric="Kappa",trControl=ctrl,tuneLength=10) # train 10 different mtry values using random search
    #===========================================================
    # apply to validation set
    pred_val = predict(m_rf,val_factors,type="prob")
    result.roc <- roc(DATval$diagn, pred_val$Control)
    aucval = result.roc$auc
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
    pred_val_cat = rep("Control",nrow(pred_val))
    pred_val_cat[pred_val$Epilepsy > 0.500] = "Epilepsy"
    confmat = makeconfmatsquare(table(pred_val_cat,LABval$diagnosis))
    result$val.confmatrix[cnt] = paste0(confmat[1,1:2],"_",confmat[2,1:2],collapse="_")
    result$val.auc[cnt] = round(aucval,digits=3)
    result$val.kappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=3)
    result$val.acc[cnt] = sum(diag(confmat)) / sum(confmat)
    #============================================================
    # testing
    pred_test = predict(m_rf, test_factors,type="prob")
    result.roc <- roc(DATtest$diagn, pred_test$Control)
    # Code to plot the ROC curve if we wanted to:
    #         x11()
    #         plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
    auctest = result.roc$auc
    result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
    pred_test_cat = rep("Control",nrow(pred_test))
    pred_test_cat[pred_test$Epilepsy > 0.500] = "Epilepsy"
    confmat = makeconfmatsquare(table(pred_test_cat,LABtest$diagnosis))
    result$test.confmatrix[cnt] = paste0(confmat[1,1:2],"_",confmat[2,1:2],collapse="_")
    result$test.auc[cnt] = round(auctest,digits=3)
    result$test.kappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=3)
    result$test.acc[cnt] = sum(diag(confmat)) / sum(confmat)
    
#   Code to plot the model error rate: randomForest:::plot.randomForest(x=m_rf$finalModel)
#     library(party)
#     x11()
#     cf = cforest(diagn ~ .,data=DATtrain[,c(fes,ncol(DATtrain))])
#     pt <- party:::prettytree(cf@ensemble[[2]], names(cf@data@get("input"))) 
#     nt <- new("BinaryTree") 
#     nt@tree <- pt 
#     nt@data <- cf@data 
#     nt@responses <- cf@responses 
#     plot(pt) 
    #https://www.analyticsvidhya.com/blog/2016/04/complete-tutorial-tree-based-modeling-scratch-in-python/
    # m_rf = train(diagn ~ .,data=DATtrain[,c(which(allvalues==valueseval[jj]),ncol(DATtrain))],method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)
    #===================================================================
    # keep track of model
    print(paste0(valueseval[jj]," validation Kappa: ",result$val.kappa[cnt],
                 " validation accuracy: ",result$val.acc[cnt],
                 " validation auc: ",result$val.auc[cnt]))
    cnt = cnt+1
  }
  if (usePCA == FALSE) {
    # select all models that are within 0.1 Kappa from the best model
    result = result[with(result,order(val.kappa)),]
    if (testparti == "wavelet") { # select only best wavelet
      bestmodels = result$model[which(result$val.kappa == max(result$val.kappa))[1]]
    } else {
      bestmodels = result$model[which(result$val.kappa > (max(result$val.kappa) - 0.2))]
    }
    if (testparti == "wavelet") {
      sel = which(modeldict$v2 %in% bestmodels == TRUE)
    } else if (testparti == "features") {
      sel = which(modeldict$v3 %in% bestmodels == TRUE)
    } else if (testparti == "waveletlevel") {
      sel = which(modeldict$v5 %in% bestmodels == TRUE)
    } else if (testparti == "aggregationtype") {
      sel = which(modeldict$v6 %in% bestmodels == TRUE)
    }
    modeldict = modeldict[sel,]
    sel = c(sel,ncol(DATtrain))
    DATtrain = DATtrain[,sel]
    DATval = DATval[,sel]
    DATtest = DATtest[,sel]
  }
  
  print(result)
  # rm(result)
}


# Now do random search to select best possible model:
# # Find out importance of parameters:
# randomvar = c(sample(1:(ncol(DATtrain)-1),10),ncol(DATtrain))
# ctrl = trainControl(method = "none") #"repeatedcv",number=3,repeats=1
# grid_rf = expand.grid(.mtry = c(2)) #2,4,8,16
# m_rf = train(diagn ~ .,data=DATtrain[,randomvar],method="rf",metric="Kappa",trControl=ctrl,tuneGrid=grid_rf)
# pred_val = predict(m_rf, DATval[,randomvar],type="prob")
# # result.roc <- roc(DATval$diagn, pred_val$Control)
# # aucval = result.roc$auc
# # result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
# # result$val.auc[cnt] = round(aucval,digits=3)
# pred_val_cat = rep("Control",nrow(pred_val))
# pred_val_cat[pred_val$Epilepsy > 0.500] = "Epilepsy"
# confmat = makeconfmatsquare(table(pred_val_cat,LABval$diagnosis))
# val.kappa = round(cohen.kappa(x=confmat)$kappa,digits=3)
# kk = modeldict[randomvar[1:(length(randomvar)-1)],]
# unique(kk$v2)
# kkk
# result$kappa = round(result$kappa,digits=3)
# result$trainingkappa  = round(result$trainingkappa,digits=3)
# write.csv(result,file="data/result_guineabissau.csv")
# fit1 = lm(val.acc ~ model,data=result)


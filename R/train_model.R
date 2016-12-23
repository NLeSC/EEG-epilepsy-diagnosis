train_model = function(DATtrain,LABtrain,DATval,LABval,modeldict,classifier="rf") { #,aggregateperid
  #only look for best possible wavelet type, but used all features and aggregationtypes
  testpart = c("wavelet") #,"features","aggregationtype") #,"waveletlevel"
  performancemetric = "Spec" #"Spec"
  for (testparti in testpart) {
    cnt = 1
    set.seed(300)
    if (testparti == "wavelet") {
      allvalues = modeldict$wvtype #v2
    } else if (testparti == "features") {
      allvalues = modeldict$feature #v3
    } else if (testparti == "waveletlevel") {
      allvalues = modeldict$wvlevel #v5
    } else if (testparti == "aggregationtype") {
      allvalues = modeldict$aggtype #v6
    }
    valueseval = unique(allvalues)
    valuerawexist = which(valueseval == "raw")
    if (length(valuerawexist) > 0) valueseval = valueseval[-valuerawexist]
    Ncases = length(valueseval)
    result = data.frame(model=rep(" ",Ncases),
                        stringsAsFactors=FALSE)
    for (jj in 1:length(valueseval)) {
      result$model[cnt] = valueseval[jj]
      
      fes = which(allvalues==valueseval[jj] | allvalues == "raw") # always add features derived from raw data
      train_factors = DATtrain[,fes]
      val_factors = DATval[,fes]
      #===========================================================
      # set training paramerters
      # ctrl = trainControl(method = "none",number=10,repeats=10,search="random")
      ctrl = caret::trainControl(method = "repeatedcv",number=10,repeats=1,search="random",
                          classProbs = TRUE, summaryFunction = caret::twoClassSummary) #, summaryFunction = twoClassSummary) # twoClassSummary is needed for Sensitivity optimization
      # ctrl = trainControl(method = "repeatedcv",number=10,repeats=3,search="grid")
      #===========================================================
      # train on training set
      set.seed(300)
      seeds <- vector(mode = "list", length = 10)#
      for(i in 1:10) seeds[[i]] <- 300
      if (classifier == "rf") {
        # random forest
        m_rf = caret::train(y=as.factor(make.names(DATtrain$diagnosis)),x=train_factors,seeds=seeds,
                     method="rf",trControl=ctrl,tuneLength=10,metric=performancemetric) # # train 10 different mtry values using random search
      }
      if (classifier == "lg") {
        m_rf = caret::train(y=as.factor(make.names(DATtrain$diagnosis)),x=train_factors,#seeds=seeds,
                     method="glm",family="binomial",trControl=ctrl,tuneLength=10,metric=performancemetric) # # train 10 different mtry values using random search
      }
      #  inspect how probabilities of classification are distributed in training set
#       pred_train = predict(m_rf,train_factors,type="prob")
#       print(summary(pred_train$X1))
#       print(quantile(pred_train$X1,probs=c(seq(0.05,0.95,by=0.05))))
#       llll
      
      #metric="Kappa"
      #metric = "Sens"
      #===========================================================
      # apply to validation set
      pred_val = predict(m_rf,val_factors,type="prob")
      #aggregate per person
      DATval_agg = DATval
      LABval_agg = LABval
      if (length(LABval$id) != length(unique(LABval$id))) {#aggregateperid == FALSE) { # aggregates estimates per person
        pred_val = data.frame(pred_val,id=LABval$id)
        pred_val = aggregate(. ~ id,data=pred_val,mean)
        DATval_agg = aggregate(. ~ id,data=DATval,mean)
        LABval_agg = aggregate(. ~ id,data=LABval,function(x){x[1]})
      }
      result.roc <- pROC::roc(DATval_agg$diagnosis, pred_val$X1) # X1 is the control group
      aucval = result.roc$auc
      result.coords <- pROC::coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
      pred_val_cat = rep("X1",nrow(pred_val))
      pred_val_cat[pred_val$X2 > 0.500] = "X2" #"Epilepsy"
      refe = make.names(LABval_agg$diagnosis)
      predi = pred_val_cat
      confmat = create_confmatrix(predi,refe)
      result$val.confmatrix[cnt] = paste0(confmat[1,1],"_",confmat[1,2],"_",confmat[2,1],"_",confmat[2,2])
      result$val.auc[cnt] = round(aucval,digits=3)
      result$val.kappa[cnt] = round(psych::cohen.kappa(x=confmat)$kappa,digits=3)
      result$val.acc[cnt] = sum(diag(confmat)) / sum(confmat)
      predi = which(names(dimnames(confmat))=="predicted")
      if (predi == 1) {  # sensitivty to detect Epilepsy
        sens = round(confmat[2,2] / (confmat[2,2]+confmat[1,2]),digits=3) 
      } else {
        sens= round(confmat[2,2] / (confmat[2,2]+confmat[2,1]),digits=3) 
      }
      result$val.sens[cnt] = sens
      #===================================================================
      print(paste0(valueseval[jj]," validation Kappa: ",result$val.kappa[cnt],
                   " validation accuracy: ",result$val.acc[cnt],
                   " validation auc: ",result$val.auc[cnt],
                   " conf matrix: ",result$val.confmatrix[cnt],
                   " validation sens: ",result$val.sens[cnt]))
      cnt = cnt+1
    }
    result = result[with(result,order(val.sens,val.kappa,val.auc,val.acc)),]
    print("Now train best model on training data and validation data combined")
    fes = which(allvalues==as.character(result$model[nrow(result)]) | allvalues == "raw")
    
    train_factors = rbind(DATtrain[,fes],DATval[,fes]) #features of train and validation combined
    din = as.factor(c(make.names(DATtrain$diagnosis),make.names(DATval$diagnosis))) #diagnosis of train and val combined
    set.seed(300)
    # train 10 different mtry values using random search
    if (classifier == "rf") {
      best_model = caret::train(y=din,x=train_factors,seeds=seeds,
                                      method="rf",metric=performancemetric,trControl=ctrl,tuneLength=10) 
    }
    if (classifier == "lg") {
      best_model = caret::train(y=din,x=train_factors,#seeds=seeds,
                                      method="glm",family="binomial",trControl=ctrl,tuneLength=10,metric=performancemetric)
    }
  }
  invisible(list(result=result,best_model=best_model,fes=fes))
}
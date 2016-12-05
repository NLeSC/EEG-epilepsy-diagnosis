train_model = function(DATtrain,LABtrain,DATval,LABval,modeldict) {
  #only look for best possible wavelet type, but used all features and aggregationtypes
  testpart = c("wavelet") #,"features","aggregationtype") #,"waveletlevel"
  for (testparti in testpart) {
    cnt = 1
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
      train_factors = DATtrain[,fes]
      val_factors = DATval[,fes]
      #===========================================================
      # set training paramerters
      # ctrl = trainControl(method = "none",number=10,repeats=10,search="random")
      ctrl = trainControl(method = "repeatedcv",number=10,repeats=1,search="random",classProbs = TRUE, summaryFunction = twoClassSummary)
      # ctrl = trainControl(method = "repeatedcv",number=10,repeats=3,search="grid")
      #===========================================================
      # train on training set

      set.seed(300)
      seeds <- vector(mode = "list", length = 10)#
      for(i in 1:10) seeds[[i]] <- 300
      m_rf = train(y=DATtrain$diagn,x=train_factors,seeds=seeds,
                   method="rf",trControl=ctrl,tuneLength=10,metric = "Sens") # # train 10 different mtry values using random search
      #metric="Kappa"
      #===========================================================
      # apply to validation set
      pred_val = predict(m_rf,val_factors,type="prob")
      result.roc <- roc(DATval$diagn, pred_val$Control)
      aucval = result.roc$auc
      result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
      pred_val_cat = rep("Control",nrow(pred_val))
      pred_val_cat[pred_val$Epilepsy > 0.500] = "Epilepsy"
      confmat = create_confmatrix(pred_val_cat,LABval$diagnosis)
      result$val.confmatrix[cnt] = paste0(confmat[1,1:2],"_",confmat[2,1:2],collapse="_")
      result$val.auc[cnt] = round(aucval,digits=3)
      result$val.kappa[cnt] = round(cohen.kappa(x=confmat)$kappa,digits=3)
      result$val.acc[cnt] = sum(diag(confmat)) / sum(confmat)
      #===================================================================
      print(paste0(valueseval[jj]," validation Kappa: ",result$val.kappa[cnt],
                   " validation accuracy: ",result$val.acc[cnt],
                   " validation auc: ",result$val.auc[cnt]))
      cnt = cnt+1
    }
    result = result[with(result,order(val.kappa,val.auc,val.acc)),]
    print("Now train best model on training data and validation data combined")
    fes = which(allvalues==as.character(result$model[nrow(result)]))
    train_factors = rbind(DATtrain[,fes],DATval[,fes]) #x
    din = as.factor(c(as.character(DATtrain$diagn),as.character(DATval$diagn))) #y
    set.seed(300)
    best_model_randomforest = train(y=din,x=train_factors,
                                    method="rf",metric="Sens",trControl=ctrl,tuneLength=10) # train 10 different mtry values using random search
  }
  invisible(list(result=result,best_model_randomforest=best_model_randomforest,fes=fes))
}
evaluatemodel = function(model,x,labels,proto_i,aggregateperid) {
  # x is evaluation data
  # labels is labels for evaluation data
  pred_test = predict(model, x,type="prob")
  labels_agg = labels
  x_agg = x
  if (aggregateperid == FALSE) {
    pred_test = data.frame(pred_test,id=labels$id)
    pred_test = stats::aggregate(. ~ id,data=pred_test,mean)
    x_agg = aggregate(. ~ id,data=x,mean) #new
    labels_agg = stats::aggregate(. ~ id,data=labels,function(x){x[1]})
  }
  # result.roc <- pROC::roc(labels_agg$diagnosis, pred_test$X1)
  result.roc <- pROC::roc(x_agg$diagnosis, pred_test$X1)
  auctest = result.roc$auc
  result.coords <- pROC::coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
  pred_test_cat = rep("X1",nrow(pred_test))
  pred_test_cat[pred_test$X2 > 0.500] = "X2"
  refe = make.names(labels_agg$diagnosis)
  predi = pred_test_cat
  # print(predi)
  # print(refe)
  confmat = create_confmatrix(predi,refe) 
  
  test.confmatrix = paste0(confmat[1,1],"_",confmat[1,2],"_",confmat[2,1],"_",confmat[2,2])
  test.auc = round(auctest,digits=3)
  test.kappa = round(psych::cohen.kappa(x=confmat)$kappa,digits=3)
  test.acc = round(sum(diag(confmat)) / sum(confmat),digits=3)
  predij = which(names(dimnames(confmat))=="predicted")
  if (length(predij) == 0) predij = 1
  if (predij == 1) {  # sensitivty to detect Epilepsy
    test.sens = round(confmat[2,2] / (confmat[2,2]+confmat[1,2]),digits=3) 
  } else {
    test.sens= round(confmat[2,2] / (confmat[2,2]+confmat[2,1]),digits=3) 
  }
  print(paste0("performance on test set: ",proto_i," acc ",test.acc," kappa ",test.kappa," auc ",test.auc," ",test.confmatrix," sens ",test.sens ))
  invisible(list(proto_i=proto_i,test.acc=test.acc,test.kappa=test.kappa,test.auc=test.auc,
            test.confmatrix=test.confmatrix,test.sens=test.sens,aggregateperid=aggregateperid))
}
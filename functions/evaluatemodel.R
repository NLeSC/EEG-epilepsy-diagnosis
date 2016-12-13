evaluatemodel = function(model,x,labels) {
  # x is evaluation data
  # labels is labels for evaluation data
  pred_test = predict(model, x,type="prob")
  labels_agg = labels
  if (perid == FALSE) {
    pred_test = data.frame(pred_test,id=labels$id)
    pred_test = aggregate(. ~ id,data=pred_test,mean)
    labels_agg = aggregate(. ~ id,data=labels,function(x){x[1]})
  }
  result.roc <- roc(labels_agg$diagnosis, pred_test$X1)
  auctest = result.roc$auc
  result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
  pred_test_cat = rep("X1",nrow(pred_test))
  pred_test_cat[pred_test$X2 > 0.500] = "X2"
  refe = make.names(labels_agg$diagnosis)
  predi = pred_test_cat
  confmat = create_confmatrix(predi,refe) 
  test.confmatrix = paste0(confmat[1,1],"_",confmat[1,2],"_",confmat[2,1],"_",confmat[2,2])
  test.auc = round(auctest,digits=3)
  test.kappa = round(cohen.kappa(x=confmat)$kappa,digits=3)
  test.acc = round(sum(diag(confmat)) / sum(confmat),digits=3)
  predi = which(names(dimnames(confmat))=="predicted")
  if (predi == 1) {  # sensitivty to detect Epilepsy
    test.sens = round(confmat[2,2] / (confmat[2,2]+confmat[1,2]),digits=3) 
  } else {
    test.sens= round(confmat[2,2] / (confmat[2,2]+confmat[2,1]),digits=3) 
  }
  print(paste0(proto_i," acc ",test.acc," kappa ",test.kappa," auc ",test.auc," ",test.confmatrix," sens ",test.sens ))
}
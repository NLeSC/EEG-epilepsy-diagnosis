reformat_DATLAB = function(DAT,LAB,aggregateperid=FALSE) {
  # turn LAB variables numberic to ease aggregation
  LAB$diagnosis_number = 2 # Epilepsy
  LAB$diagnosis_number[which(LAB$diagnosis == "Control")] = 1
  LAB$protocol_number = 2 #eyesclosed
  LAB$protocol_number[which(LAB$protocol == "eyesopen")] = 1
  LAB = LAB[,-c(which(names(LAB)=="fnames" | names(LAB)=="protocol" | names(LAB)=="diagnosis"))]
  names(LAB)[which(names(LAB)=="diagnosis_number")] = "diagnosis"
  names(LAB)[which(names(LAB)=="protocol_number")] = "protocol"
  LAB$fnames = rownames(LAB)
  
  if (aggregateperid == TRUE) {
    # aggregate DAT multiple windows and epochs per unique identifier id
    DATLAB = merge(DAT,LAB,by="fnames")
    DATLAB = DATLAB[,-c(1,(ncol(DATLAB)-3):(ncol(DATLAB)-2))]
    DATLAB$id = as.factor(DATLAB$id)
    DAT = aggregate(. ~ id + diagnosis + protocol,data=DATLAB,mean)
    # tidy up DAT
    DAT = DAT[,-which(names(DAT)=="window" | names(DAT)=="fnames")]
    movecolumn = which(names(DAT)=="id" | names(DAT)=="diagnosis" | names(DAT) == "protocol")
    DAT = cbind(DAT,DAT[,movecolumn])
    DAT = DAT[,-movecolumn] # delete double columns
    DAT$id = as.numeric(as.character(DAT$id))
    # aggregate LAB
    LAB = aggregate(. ~ id + diagnosis + protocol,data=LAB,function(x){x[1]})
    
  } else {
    DATLAB = merge(DAT,LAB,by="fnames")
    DAT = DATLAB[,-which(names(DATLAB)=="window" | names(DATLAB) == "dur" | names(DATLAB) == "epoch" | names(DATLAB) == "fnames")]
    nDB = names(DATLAB)
    LAB = DATLAB[,nDB[which(nDB %in% names(LAB) == TRUE)]]
  }
  invisible(list(DAT=DAT,LAB=LAB))
}
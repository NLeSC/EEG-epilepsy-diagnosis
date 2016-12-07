split_data = function(LAB,DAT,logfile = "log_guinneabissau.csv",proto_i= "eyesopen",split=c(20,20),
                      uselog = TRUE,logdur=4,useallepoch=FALSE) {
  
  split_data = function(LAB,DAT,proto_i,split=c(20,20)) {
    ids = c()
    if (useallepoch == TRUE) {
      ep = LAB$epoch >= 1
    } else {
      ep = LAB$epoch == 1
    }
    
    for (set in 1:2) { # first select test and validation set based on split as provided in the arguments
      for (diag_i in c("Control","Epilepsy")) {
        x = LAB[which(ep & LAB$protocol == proto_i & LAB$diagnosis == diag_i & 
                        (LAB$id %in% ids) == FALSE),]$id
        if (length(x) == 0) {
          # investigate what info is missing
          if (length(which(LAB$protocol == proto_i)) == 0) print("requested protocol not in dataset")
          if (length(which(ep == TRUE)) == 0) print("no epochs in dataset")
          if (length(which(LAB$diagnosis == diag_i)) == 0) print("patient group not in dataset")
          if (length(which((LAB$id %in% ids) == FALSE)) == 0) print("no id left")
        }
        set.seed(300)
        ids = c(ids,sample(x,size=round(split[set]/2)))
      }
      if (set == 1) { # validation
        LABval = LAB[which(ep & LAB$protocol == proto_i & (LAB$id %in% ids)==TRUE),]
        DATval = DAT[which(DAT$fnames %in% rownames(LABval)  == TRUE),]
        ids_val = ids
      } else { # test
        ids = ids[which(ids %in% ids_val == FALSE)]
        LABtest = LAB[which(ep & LAB$protocol == proto_i & (LAB$id %in% ids)==TRUE),] #ep == 1
        DATtest = DAT[which(DAT$fnames %in% rownames(LABtest)  == TRUE),]
        ids_test = ids
      }
    }
    # all remaining data will go to the training set
    # ids = ids[which(ids %in% ids_test == FALSE)]
    LABtrain = LAB[which(ep & LAB$protocol == proto_i & (LAB$id %in% ids_val)==FALSE
                         & (LAB$id %in% ids_test)==FALSE),]
    DATtrain= DAT[which(DAT$fnames %in% rownames(LABtrain)  == TRUE),]
    invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain))
  }
  split_data_bypython = function(logfile,proto_i,logdur) {
    # this log comes from the python notebook, this is where the data is split up!
    LOG = read.csv(logfile,stringsAsFactors = FALSE)
    if (length(which(LOG$protocol == proto_i)) == 0) {
      warning("value of argument proto_i was not recognized in the log file")
    }
    relevantfilenames = function(LOG,setname) {
      fn = LOG[which(LOG$set == setname & LOG$dur == logdur & LOG$protocol == proto_i),]$filename
      return(fn)
    }
    train_fnames = relevantfilenames(LOG,setname="train")
    traini = which((rownames(LAB) %in% train_fnames) == TRUE)
    LABtrain = LAB[traini,]
    DATtrain = DAT[traini,]
    # test set
    test_fnames = relevantfilenames(LOG,setname="test")
    testi = which((rownames(LAB) %in% test_fnames) == TRUE)
    LABtest = LAB[testi,]
    DATtest = DAT[testi,]
    # validation set
    val_fnames = relevantfilenames(LOG,setname="valid")
    vali = which((rownames(LAB) %in% val_fnames) == TRUE)
    LABval = LAB[vali,]
    DATval = DAT[vali,]
    invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain))
  }
  # make sure that both LAB and DAT have matching row order
  getid = function(x) {
    tmp = unlist(strsplit(x,"_"))[2]
    tmp2 = unlist(strsplit(tmp,"id"))[2]
    return(as.numeric(tmp2))
  }
  #-----------------------------------------------------
  DAT$id = as.numeric(sapply(DAT$fnames,getid))
  DAT = DAT[order(DAT$id),]
  if (length(which((rownames(LAB) == unique(DAT$fnames)) == FALSE)) > 0) print("Error: data order does not match")
  
  
  if (uselog == TRUE) {
    P = split_data_bypython(logfile,proto_i,logdur,split,logdur)
    LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  } else {
    P = split_data(LAB,DAT,proto_i,split=split) #"eyesopen" #,"eyesclosed"
    LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  }
  DATtest$diagn = as.factor(LABtest$diagnosis)
  DATtrain$diagn = as.factor(LABtrain$diagnosis)
  DATval$diagn = as.factor(LABval$diagnosis)
  invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain)) 
}
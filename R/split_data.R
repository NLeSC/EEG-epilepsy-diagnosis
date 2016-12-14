split_data = function(LAB,DAT,logfile = "log_guinneabissau.csv",proto_i= 1,split=c(20,20),
                      uselog = TRUE,logdur=4) { #,sampleidentifier="id"
  # protocol 1 = eyes open
  # protocol 2 = eyes closed
  # diagnosis 1 = Control
  # diagnosis 2 = Epilepsy
  
  split_data = function(LAB,DAT,proto_i,split=c(20,20)) {
    ids = c() # the identifiers, could be id or fnames
    for (set in 1:2) { # first select test and validation set based on split as provided in the arguments
      for (diag_i in c(1,2)) { #c("Control","Epilepsy")) {
        x = unique(LAB[which(LAB$protocol == proto_i & LAB$diagnosis == diag_i & (LAB$id %in% ids) == FALSE),]$id)
        if (length(x) == 0) {
          # investigate what info is missing
          if (length(which(LAB$protocol == proto_i)) == 0) print("requested protocol not in dataset")
          # if (length(which(ep == TRUE)) == 0) print("no epochs in dataset")
          if (length(which(LAB$diagnosis == diag_i)) == 0) print("patient group not in dataset")
          if (length(which((LAB$id %in% ids) == FALSE)) == 0) print("no id left")
        }
        set.seed(300)
        ids = c(ids,sample(x,size=round(split[set]/2)))
      }
      if (set == 1) { # validation
        LABval = LAB[which(LAB$protocol == proto_i & (LAB$id %in% ids) ==TRUE),]
        DATval = DAT[which(DAT$protocol == proto_i & (DAT$id %in% ids) == TRUE),]
        ids_val = ids
      } else if (set == 2) { # test
        ids = ids[which(ids %in% ids_val == FALSE)]
        LABtest = LAB[which(LAB$protocol == proto_i & (LAB$id %in% ids)==TRUE),] #ep == 1
        DATtest = DAT[which(DAT$protocol == proto_i & (DAT$id %in% ids)  == TRUE),]
        ids_test = ids
      }
    }
    # all remaining data will go to the training set
    LABtrain = LAB[which(LAB$protocol == proto_i & ((LAB$id %in% ids_val)==FALSE & (LAB$id %in% ids_test)==FALSE)),]
    DATtrain= DAT[which(DAT$protocol == proto_i & (DAT$id %in% LABtrain$id)  == TRUE),]
    
    invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain))
  }
  #   split_data_bypython = function(logfile,proto_i,logdur) {
  #     # this log comes from the python notebook, this is where the data is split up!
  #     LOG = read.csv(logfile,stringsAsFactors = FALSE)
  #     if (length(which(LOG$protocol == proto_i)) == 0) {
  #       warning("value of argument proto_i was not recognized in the log file")
  #     }
  #     relevantfilenames = function(LOG,setname) {
  #       fn = LOG[which(LOG$set == setname & LOG$dur == logdur & LOG$protocol == proto_i),]$filename
  #       return(fn)
  #     }
  #     train_fnames = relevantfilenames(LOG,setname="train")
  #     traini = which((rownames(LAB) %in% train_fnames) == TRUE)
  #     LABtrain = LAB[traini,]
  #     DATtrain = DAT[traini,]
  #     # test set
  #     test_fnames = relevantfilenames(LOG,setname="test")
  #     testi = which((rownames(LAB) %in% test_fnames) == TRUE)
  #     LABtest = LAB[testi,]
  #     DATtest = DAT[testi,]
  #     # validation set
  #     val_fnames = relevantfilenames(LOG,setname="valid")
  #     vali = which((rownames(LAB) %in% val_fnames) == TRUE)
  #     LABval = LAB[vali,]
  #     DATval = DAT[vali,]
  #     invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain))
  #   }
  # make sure that both LAB and DAT have matching row order
#   getid = function(x) {
#     tmp = unlist(strsplit(x,"_"))[2]
#     tmp2 = unlist(strsplit(tmp,"id"))[2]
#     return(as.numeric(tmp2))
#   }
  #-----------------------------------------------------
  # DAT$id = as.numeric(sapply(DAT$fnames,getid))
  # DAT = DAT[order(DAT$id),]
  # if (length(which((rownames(LAB) == unique(DAT$id)) == FALSE)) > 0) print("Error: data order does not match")
  if (length(which(LAB$id %in% DAT$id == TRUE)) != length(DAT$id)) print("Error: data order does not match")
  
  
  if (uselog == TRUE) {
    #     P = split_data_bypython(logfile,proto_i,logdur,split,logdur)
    #     LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  } else {
    P = split_data(LAB,DAT,proto_i,split=split) #"eyesopen" #,"eyesclosed"
    LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  }
  invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain)) 
}
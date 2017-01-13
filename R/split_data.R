split_data = function(LAB,DAT,proto_i= 1,split=c(20,20),seed) {
  # protocol 1 = eyes open
  # protocol 2 = eyes closed
  # diagnosis 1 = Control
  # diagnosis 2 = Epilepsy
  split_data = function(LAB,DAT,proto_i,split=c(20,20),seed) {
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
        set.seed(seed)
        ids = c(ids,sample(x,size=round(split[set]/2)))
      }
      if (set == 1) { # validation
        LABval = LAB[which(LAB$protocol == proto_i & (LAB$id %in% ids) == TRUE),]
        DATval = DAT[which(DAT$protocol == proto_i & (DAT$id %in% ids) == TRUE),]
        ids_val = ids
      } else if (set == 2) { # test
        ids = ids[which(ids %in% ids_val == FALSE)]
        LABtest = LAB[which(LAB$protocol == proto_i & (LAB$id %in% ids) == TRUE),] #ep == 1
        DATtest = DAT[which(DAT$protocol == proto_i & (DAT$id %in% ids)  == TRUE),]
        ids_test = ids
      }
    }
    # all remaining data will go to the training set
    LABtrain = LAB[which(LAB$protocol == proto_i & ((LAB$id %in% ids_val)==FALSE & (LAB$id %in% ids_test)==FALSE)),]
    DATtrain= DAT[which(DAT$protocol == proto_i & (DAT$id %in% LABtrain$id)  == TRUE),]
    
    invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain))
  }
  #-----------------------------------------------------
  if (length(which(LAB$id %in% DAT$id == TRUE)) != length(DAT$id)) print("Error: data order does not match")
  P = split_data(LAB,DAT,proto_i,split=split,seed=seed) #"eyesopen" #,"eyesclosed"
  LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  invisible(list(LABval=LABval,LABtest=LABtest,LABtrain=LABtrain,DATval=DATval,DATtest=DATtest,DATtrain=DATtrain)) 
}
create_featuredict = function(DAT) {
  namesofplifeatures = c("meanpli","sdpli","ninetyppli","densitypli","maxdegreepli",
                    "maxstrengthpli","mpli","diameterpli","leafnumberpli","maxbcpli",
                         "eccpli","radiuspli","Thpli","kappapli")
    
  getfeaturetype = function(x) {
    tt = unlist(strsplit(x,"[.]"))
    if (length(tt) >= 4) {
      tt = tt[2]
    } else {
      tt = tt[1]
    }
    return(tt)
  }
  getwavelettype = function(x) {
    tt = unlist(strsplit(x,"[.]"))
    if (length(tt) >= 4) {
      tt = tt[length(tt)-1]
    } else {
      if (length(which(tt[1] %in% namesofplifeatures == TRUE)) > 0) {
        tt = tt[2]
      } else {
        tt = "na"
      }
    }
    return(tt)
  }
  getwaveletlevel= function(x) {
    tt = unlist(strsplit(x,"[.]"))
    if (length(tt) >= 4) {
      tt = tt[1]
    } else {
      if (length(which(tt[1] %in% namesofplifeatures == TRUE)) > 0) {
        tt = tt[3]
        if (tt == "notapplicable") tt = 0
      } else {
        tt = "na"
      }
    }
    return(tt)
  }
  getaggtype= function(x) {
    tt = unlist(strsplit(x,"[.]"))
    if (length(tt) >= 4) {
      tt = tt[4]
    } else {
      tt = "na"
    }
    return(tt)
  }
  #---------------------------
  DAT = DAT[,-which(names(DAT) =="id" | names(DAT) =="diagnosis" | names(DAT) =="protocol")] #,ncol(DAT)-1
  varnames = names(DAT)
  featuredict = data.frame(wvtype = sapply(varnames,FUN = getwavelettype), #v2
                         feature = sapply(varnames,FUN = getfeaturetype), #v3
                         wvlevel = sapply(varnames,FUN = getwaveletlevel), #v5
                         aggtype = sapply(varnames,FUN = getaggtype),stringsAsFactors=FALSE) #v6
  return(featuredict)
}
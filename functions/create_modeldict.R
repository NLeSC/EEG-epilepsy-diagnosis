create_modeldict = function(DAT) {
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
  #---------------------------
  DAT = DAT[,-c(ncol(DAT))] #,ncol(DAT)-1
  varnames = names(DAT)
  modeldict = data.frame(v2 = sapply(varnames,FUN = getwavelettype),
                         v3 = sapply(varnames,FUN = getfeaturetype),
                         v5 = sapply(varnames,FUN = getwaveletlevel),
                         v6 = sapply(varnames,FUN = getaggtype),stringsAsFactors=FALSE)
  return(modeldict)
}
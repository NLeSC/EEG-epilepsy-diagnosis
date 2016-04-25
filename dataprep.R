rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
datadir = "data/eeg"
labels = read.csv("data/labels.csv")
require(reshape); require(entropy)
require(randomForest); require(signal); require(moments)
# ; require(pracma); require(fractal)

cn = c('Fp2', 'Fp1', 'F8', 'F4', 'Fz', 'F3', 'F7', 
       'A2','A1','T8','T7','C4', 'Cz', 'C3',  
       'P8', 'P4', 'Pz', 'P3', 'P7','O2', 'O1') # channel names
fn = c("mean","sd","entropy","en.entropy","max","min","skewness","median","domfreq","maxpow") # feature names ,"lyapunov"

#--------------------------------------------------------
# define generic function for feature extraction
# Note: approx_entropy from pracma package is very slow
getfeatures = function(x,fn) {
  # x: dataframe with multiple signals in columns over which the feature is calculated
  # fn: feature names that need to be calculated
  en.entropy = function(x) { # entropy
    b = discretize(x,numBins=100)
    en.entropy = entropy::entropy(b)
  }
  domfreq = function(x) {
    pp = spectrum(x,plot=FALSE)
    domfreq = pp$freq[which.max(pp$spec)]
  }
  maxpow = function(x) {
    pp = spectrum(x,plot=FALSE)
    maxpow = max(pp$spec^2)
  }
  t = c()
  if (length(which(fn=="mean")) >0) t = rbind(t,sapply(x,mean))
  if (length(which(fn=="sd")) >0) t = rbind(t,sapply(x,sd))
  if (length(which(fn=="entropy")) >0) t = rbind(t,sapply(x,entropy))
  if (length(which(fn=="en.entropy")) >0) t = rbind(t,sapply(x,en.entropy))
  if (length(which(fn=="max")) >0) t = rbind(t,sapply(x,max))
  if (length(which(fn=="min")) >0) t = rbind(t,sapply(x,min))
  if (length(which(fn=="skewness")) >0) t = rbind(t,sapply(x,moments::skewness))
  if (length(which(fn=="median")) >0) t = rbind(t,sapply(x,median))
  # if (length(which(fn=="lyapunov")) >0) t = rbind(t,sapply(x,lyapunov))
  if (length(which(fn=="domfreq")) >0) t = rbind(t,sapply(x,domfreq))
  if (length(which(fn=="maxpow")) >0) t = rbind(t,sapply(x,maxpow))
  A = as.data.frame(t)
}
#--------------------------------------------------------
# extract patient identifiers from filename
fnames = list.files(datadir,include.dirs=FALSE,full.names = FALSE)
myfun = function(x) {
  myfun = as.numeric(unlist(strsplit(x,"_ep"))[1])
}
idnames = sapply(fnames,myfun)
uid = unique(idnames)
#--------------------------------------------------------
# Load data and extract features
files = list.files(datadir,include.dirs=TRUE,full.names = TRUE)
rowsd = rowmed = S = c()
cnt = 0
for (i in 1:length(uid)) { #unique id numbers
  if (cnt == 10) {
    prog = round((i/length(uid)) * 1000)/ 10
    print(paste0(prog," %"))
    cnt = 0
  }
  cnt = cnt + 1
  ind = which(idnames == i) #which epochs belong to this person
  if (length(ind) > 0) {
    data = read.table(files[ind[1]],sep="") # only look at first 4 second epoch for now
    if (nrow(data) >= 4096) { # only use files with at least 4096 data points
      data = data[1:4096,]
      sc = t(data)
      #--------------------------------------------
      # extract wavelets:
      require(wavelets)
      wtdata = NULL
      for (i in 1:nrow(sc)) {
        a <- as.matrix(sc[i,]) #t(sc[i,])
        wt <- dwt(a, filter="d4", boundary="periodic")
        wtdata <- rbind(wtdata, unlist(c(wt@W,wt@V[[wt@level]])))
      }
      wtdata <- as.data.frame(wtdata)# all the 12 wavelets that are possible in a 4607 signal
      
      bands = c(rep(1,2048),rep(2,1024),rep(3,512),rep(4,256),rep(5,128),rep(6,64),rep(7,32),
                rep(8,16),rep(9,8),rep(10,4),rep(11,2),rep(12,1),rep(0,1)) # see also names(wtdata)
      # extract features from wavelets:
      for (j in 1:9) {
        df = as.data.frame(t(wtdata[,which(bands==j)]))
        A = getfeatures(df,fn)
        colnames(A) = cn
        # also calculate median feature values
        tmp = c(as.numeric(apply(A,1,median)),idnames[ind[1]],j)
        tmp = as.data.frame(t(tmp))
        # also calculate sd feature values
        tmp2 = c(as.numeric(apply(A,1,sd)),idnames[ind[1]],j)
        tmp2 = as.data.frame(t(tmp2))
        # add labels and merge
        colnames(tmp) = colnames(tmp2) = c(fn,"id","band")
        rowmed = rbind(rowmed,tmp)
        rowsd = rbind(rowsd,tmp2)
        A$method = paste0(fn,j)
        A$id = idnames[ind[1]]
        S = rbind(S,A)
      }
      #--------------------------------------------
      # extract features from the raw signals:
      A = getfeatures(data,fn)
      colnames(A) <- cn
      
      # also calculate median feature values
      tmp = c(as.numeric(apply(A,1,median)),idnames[ind[1]],0)
      tmp = as.data.frame(t(tmp))
      # also calculate median feature values
      tmp2 = c(as.numeric(apply(A,1,sd)),idnames[ind[1]],0)
      tmp2 = as.data.frame(t(tmp2))
      # add labels and merge
      colnames(tmp) = colnames(tmp2) = c(fn,"id","band")
      rowmed = rbind(rowmed,tmp)
      rowsd = rbind(rowsd,tmp2)
      A$method = paste0(fn," raw")
      A$id = idnames[ind[1]]
      S = rbind(S,A)
      rm(A)
    }
  }
}
# summary of features across channels
rowsd2 = reshape(data=rowsd,timevar="band",idvar="id",direction="wide")
rowsd3 = rowsd2[which(rowsd2$id %in% labels$DB_number == TRUE),]
DATsd = rowsd3[,-which(colnames(rowsd3)=="id")]
rowmed2 = reshape(data=rowmed,timevar="band",idvar="id",direction="wide")
rowmed3 = rowmed2[which(rowmed2$id %in% labels$DB_number == TRUE),]
DATmed = rowmed3[,-which(colnames(rowmed3)=="id")]
# features per channels
W = reshape(data=S,timevar="method",idvar="id",direction="wide")
W2 = W[which(W$id %in% labels$DB_number == TRUE),]
DAT = W2[,-which(colnames(W2)=="id")]
LAB = labels[which(labels$DB_number %in% W2$id == TRUE),]
save(DAT,DATmed,DATsd,LAB,labels,file="features.RData")

# approx.entropy = function(ts) { # aproximate entropy
#   approx.entropy = pracma::approx_entropy(ts, edim = 2, r = 0.2*sd(ts), elag = 1)
# }
# bf.fil = function(x) { # band-pass filter
#   sf = 1000
#   lb = 0.2
#   hb = 50
#   Wc = matrix(0,2,1) 
#   Wc[1,1] = lb / (sf/2) 
#   Wc[2,1] = hb / (sf/2) 
#   n = 4
#   bf = signal::butter(n,Wc,type=c("pass")) 
#   bf.fil = signal::filter(bf,x) 
# }

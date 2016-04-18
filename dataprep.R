rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
datadir = "data/eeg"
labels = read.csv("data/labels.csv")
require(reshape); require(entropy)
require(randomForest); require(signal); require(pracma)

# en.entropy = function(x) { # entropy
#   b = discretize(x,numBins=100)
#   en.entropy = entropy::entropy(b)
# }
#========================================
# extract patient identifiers from filename
fnames = list.files(datadir,include.dirs=FALSE,full.names = FALSE)
myfun = function(x) {
  myfun = as.numeric(unlist(strsplit(x,"_ep"))[1])
}
idnames = sapply(fnames,myfun)
uid = unique(idnames)
#========================================
# Load data and extract features
files = list.files(datadir,include.dirs=TRUE,full.names = TRUE)
S = c()
for (i in 1:length(uid)) { #unique id numbers
  ind = which(idnames == i) #which epochs belong to this person
  if (length(ind) > 0) {
    data = read.table(files[ind[1]],sep="") # only look at first 4 second epoch for now
    if (nrow(data) >= 4096) { # only use files with at least 4096 data points
      data = data[1:4096,]
      sc = t(data)
      require(wavelets)
      wtdata = NULL
      for (i in 1:nrow(sc)) {
        a <- as.matrix(sc[i,]) #t(sc[i,])
        wt <- dwt(a, filter="haar", boundary="periodic")
        wtdata <- rbind(wtdata, unlist(c(wt@W,wt@V[[wt@level]])))
      }
      wtdata <- as.data.frame(wtdata)# all the 12 wavelets that are possible in a 4607 signal
      bands = c(rep(1,2048),rep(2,1024),rep(3,512),rep(4,256),rep(5,128),rep(6,64),rep(7,32),
                rep(8,16),rep(9,8),rep(10,4),rep(11,2),rep(12,1),rep(0,1)) # see also names(wtdata)
      for (j in 1:11) {
        df = as.data.frame(t(wtdata[,which(bands==j)]))
        A = as.data.frame(rbind(sapply(df,mean),sapply(df,sd),sapply(df,entropy))) #,sapply(df,en.entropy),sapply(df,approx_entropy)
        colnames(A) <- c( 'Fp2', 'Fp1', 'F8', 'F4', 'Fz', 'F3', 'F7', 
                          'A2','A1','T8','T7','C4', 'Cz', 'C3',  
                          'P8', 'P4', 'Pz', 'P3', 'P7','O2', 'O1' )
        A$method = c(paste0("mean",j),paste0("sd",j),paste0("entropy",j)) #"median","entropy",
        A$id = idnames[ind[1]]
        S = rbind(S,A)
      }
      rm(A)
    }
  }
}
W = reshape(data=S,timevar="method",idvar="id",direction="wide")
W2 = W[which(W$id %in% labels$DB_number == TRUE),]
DAT = W2[,-which(colnames(W2)=="id")]
LAB = labels[which(labels$DB_number %in% W2$id == TRUE),]
save(DAT,LAB,labels,file="features.RData")

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

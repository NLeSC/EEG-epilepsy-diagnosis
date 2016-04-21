rm(list=ls())
graphics.off()

setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
datadir = "data/eeg"

labels = read.csv("data/labels.csv")
require(reshape); require(entropy)
require(randomForest); require(signal); require(pracma)

#========================================
# extract patient identifiers from filename
fnames = list.files(datadir,include.dirs=FALSE,full.names = FALSE)
myfun = function(x) {
  myfun = as.numeric(unlist(strsplit(x,"_ep"))[1])
}
idnames = sapply(fnames,myfun)
uid = unique(idnames)
#========================================
# Load data and visualise
files = list.files(datadir,include.dirs=TRUE,full.names = TRUE)
S = c()
yesep = noep = 0
require(wavelets)
for (i in 1:length(uid)) { #unique id numbers
  ind = which(idnames == i) #which epochs belong to this person
  if (length(ind) > 0) {
    data = read.table(files[ind[1]],sep="") # only look at first 4 second epoch for now
    if (nrow(data) >= 4096) { # only use files with at least 4096 data points
      data = data[1:4096,]
      
      #       #----------------------
#       sc = t(data)
#       # extract wavelets:
#       wtdata = NULL
#       for (i in 1:nrow(sc)) {
#         a <- as.matrix(sc[i,]) #t(sc[i,])
#         wt <- dwt(a, filter="haar", boundary="periodic")
#         wtdata <- rbind(wtdata, unlist(c(wt@W,wt@V[[wt@level]])))
#       }
#       wtdata <- as.data.frame(wtdata)# all the 12 wavelets that are possible in a 4607 signal
#       bands = c(rep(1,2048),rep(2,1024),rep(3,512),rep(4,256),rep(5,128),rep(6,64),rep(7,32),
#                 rep(8,16),rep(9,8),rep(10,4),rep(11,2),rep(12,1),rep(0,1)) # see also names(wtdata)
#       # extract features from wavelets:
#       # for (j in 1:11) {
#       df = as.data.frame(t(wtdata[,which(bands==7)]))
#       data = df
      #----------------
      t2 = 4000
      if (t2 > nrow(data)) t2 = nrow(data)
      yrange = c(-80,80) #c(0,10^5)
      xrange =c(0,100) #c(-2,2)
      id = idnames[ind[1]]
      info = labels[which(labels$DB_number == id),] #1 no epilepsy, 2 epilepsy
      myfun = function(x) {
        tmp  =x
        myfun = tmp/sd(tmp)
#         tmp  =diff(x)
#         myfun = tmp/sd(tmp)
      }
      if (info$definitive_diagnosis == 1 & noep < 4 | info$definitive_diagnosis == 2 & yesep < 4) {
        if (info$definitive_diagnosis == 1) {
          coll =  "black"
        } else {
          coll = "blue"
        }
        x11()
        par(mfrow=c(7,3),oma=c(1,1,1,1),mar=c(1,1,1,1))
        psd = c()
        for (j in 1:21) {
          # plot(myfun(x=data[1:t2,j]),type="l",ylim=yrange,col=coll)
          plot(data[1:t2,j],type="l",ylim=yrange,col=coll)
#           sp = spectrum(data[1:t2,j],plot=TRUE,col=coll)
#           psd = rbind(psd,sp$spec)
        }
        # plot(myfun(x=rowMeans(data[1:t2,])),type="l",ylim=yrange,col=coll)
        # psdsum = colMeans(psd)
        # plot(psdsum,type="l",ylim=yrange,col=coll) #log(
        if (noep < 4 & info$definitive_diagnosis == 1) noep = noep + 1
        if (yesep < 4 & info$definitive_diagnosis == 2) yesep = yesep + 1
      }
      if (yesep == 3 & noep == 3) jjjj
    }
  }
}
# W = reshape(data=S,timevar="method",idvar="id",direction="wide")
# W2 = W[which(W$id %in% labels$DB_number == TRUE),]
# DAT = W2[,-which(colnames(W2)=="id")]
# LAB = labels[which(labels$DB_number %in% W2$id == TRUE),]

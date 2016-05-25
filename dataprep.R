rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
source("getfeatures.R")
datadir = "data/eeg"
labels = read.csv("data/labels.csv")
require(reshape)
cn = c('Fp2', 'Fp1', 'F8', 'F4', 'Fz', 'F3', 'F7', 
       'A2','A1','T8','T7','C4', 'Cz', 'C3',  
       'P8', 'P4', 'Pz', 'P3', 'P7','O2', 'O1') # channel names
# featurenames to be extracted in a vector:
# fn = c("mean","sd","entropy","en.entropy","max","min","skewness",
#        "median","domfreq","maxpow","zerocross","RMS") # feature names ,"lyapunov"
fn = c("sd","entropy","max","min")

lf.fil = function(x) { # low-pass filter
  sf = 1000
  lb = 60
  Wc = lb / (sf/2) 
  n = 4
  lf = signal::butter(n,Wc,type=c("low")) 
  lf.fil = signal::filter(lf,x) 
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
S = c()
cnt = 0
# filtertypes =  c(paste0("d",seq(2,20,by=2)), # Daubechies
#                  paste0("la",seq(8,20,by=2)), #Least Aymetric
#                  paste0("bl",c(14,18,20)), #Best localized
#                  paste0("c",seq(6,30,by=6))) # Coiflet
filtertypes =  c(paste0("d",seq(2,20,by=2)), # Daubechies
                 paste0("la",seq(8,20,by=2)), #Least Aymetric
                 paste0("bl",c(14,18,20)), #Best localized
                 paste0("c",seq(6,30,by=6))) # Coiflet

require(wavelets)
# Wavelets (?):
#1: 1024-2048 Hz;  #2: 512-1024 Hz; #3: 256-512 Hz
#4: 128-256 Hz;    #5: 64-128 Hz;   #6: 32-64 Hz
#7: 16-32 Hz beta; #8: 8-16 Hz alpha;    #9: 4-8 Hz theta
#10: 2-4 Hz delta; #11: 1-2 Hz delta;    #12: 0.5-1 Hz delta
print(paste("filtertypes: ",paste(filtertypes,collapse=" ")))
print(paste("N persons: ",paste(length(uid))))
for (i in 1:length(uid)) { #unique id numbers
  if (cnt == 2) {
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
      wtdata = NULL
      G = c()
      for (firi in filtertypes) { #filter types
        mymra = function(x){
          out = mra(x,filter=firi, boundary="periodic",n.levels=12)
          return(unlist(c(out@D)))
        }
        wtdata = t(apply(sc,1,mymra))
        wtdata <- as.data.frame(wtdata)
        bands = c(rep(1:12,each=4096))  
        temp = as.data.frame(t(wtdata))
        temp$bands = bands
        for (featuresi in fn) { # features
          A = getfeatures(x=temp,fns=featuresi)
          colnames(A) = cn
          tmp = as.numeric(apply(A,1,min)) # minimum features value across channels
          tmp = as.data.frame(t(tmp))
          tmp2 = as.numeric(apply(A,1,max)) # maximum features value across channels
          tmp2 = as.data.frame(t(tmp2))
          tmp3 = cbind(tmp,tmp2)
          colnames(tmp3) =c(paste0(1:12,".",featuresi,".",firi,".min"),paste0(1:12,".",featuresi,".",firi,".max")) # add labels and merge
          if (length(G) == 0) {
            G = tmp3
          } else {
            G = cbind(G,tmp3)  
          }
        }
        rm(tmp3)
      }
      G$id = idnames[ind[1]]
      S = rbind(S,G)
      rm(G)
    }
  }
  # print("summary characteristics of features across channels")
  # S2 = reshape(data=S,timevar="band",idvar="id",direction="wide")  
}
S2 = S
S3 = S2[which(S2$id %in% labels$DB_number == TRUE),]
DAT = S3[,-which(colnames(S3)=="id")]
LAB = labels[which(labels$DB_number %in% S2$id == TRUE),]
save(DAT,LAB,labels,file=paste0("data/features.RData")) #DAT,
# approx.entropy = function(ts) { # aproximate entropy
#   approx.entropy = pracma::approx_entropy(ts, edim = 2, r = 0.2*sd(ts), elag = 1)
# }
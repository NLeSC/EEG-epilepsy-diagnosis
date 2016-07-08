rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
source("getfeatures.R")
datadir = "/media/windows-share/EEGs_Guinea-Bissau_cleaned"
# require(reshape)
library(wavelets)
#==========================================================
# featurenames and wavelet filter types to be extracted:
# fn = c("mean","sd","entropy","en.entropy","max","min","skewness",
#        "median","domfreq","maxpow","zerocross","RMS") # feature names ,"lyapunov"
fn = c("sd","entropy","max","min") #
# filtertypes =  c(paste0("d",seq(2,20,by=2)), # Daubechies
#                  paste0("la",seq(8,20,by=2)), #Least Aymetric
#                  paste0("bl",c(14,18,20)), #Best localized
#                  paste0("c",seq(6,30,by=6))) # Coiflet
filtertypes =  c(paste0("d",seq(2,20,by=2)), # Daubechies
                 paste0("la",seq(8,20,by=2)), #Least Aymetric
                 paste0("bl",c(14,18,20)), #Best localized
                 paste0("c",seq(6,30,by=6))) # Coiflet
#==========================================================
# other parameters
sf = 128
n.levels = 7
siglen = 30*sf
# filenames
files = list.files(datadir,include.dirs=TRUE,full.names = TRUE)
files_short = list.files(datadir,include.dirs=FALSE,full.names = FALSE)
# extract metadata from filenames
getid = function(x) {
  tmp = unlist(strsplit(x,"_"))[2]
  tmp2 = unlist(strsplit(tmp,"id"))[2]
  return(as.numeric(tmp2))
}
getprotocol = function(x) return(unlist(strsplit(x,"_"))[1])
getdiagnosis = function(x) {
  tmp = unlist(strsplit(x,"_"))[4]
  tmp2 = unlist(strsplit(tmp,"gro"))[2]
  return(unlist(strsplit(tmp2,"[.]"))[1])
}
getdur = function(x) {
  tmp = unlist(strsplit(x,"_"))[3]
  tmp2 = unlist(strsplit(tmp,"dur"))[2]
  return(as.numeric(tmp2))
}
metadata = data.frame(id = as.numeric(sapply(files_short,getid)),
                      protocol = as.character(sapply(files_short,getprotocol)),
                      diagnosis =  as.character(sapply(files_short,getdiagnosis)),
                      dur = as.numeric(sapply(files_short,getdur)),row.names=files_short,stringsAsFactors = FALSE)
metadata = metadata[order(metadata$id),]
# initialize other parameters:
S = c() #object in which features will be collected
cnt = 0 #counter for showing process in console
print(paste("filtertypes: ",paste(filtertypes,collapse=" ")))
print(paste("N files: ",paste(length(files_short))))

count = rep(0,4)

for (i in 1:length(files_short)) { #unique id numbers #length(files_short)
  if (cnt == 4) { # print progress of processing
    prog = round((i/length(files_short)) * 1000)/ 10
    print(paste0(prog," %"))
    cnt = 0
  }
  cnt = cnt + 1
  data = read.csv(files[i]) # load files
  if (nrow(data) >= siglen) { # only use first 4 seconds
    data = data[1:siglen,2:15] # we are only interested in these columns
    
    #----------------
    t2 = siglen
    if (t2 > nrow(data)) t2 = nrow(data)
    yrange = c(-2,2) #c(0,10^5)
    xrange =c(0,siglen) #c(-2,2)
    myfun = function(x) {
      tmp  =x - median(x)
      myfun = tmp/sd(tmp)
    }
    info = metadata[which(rownames(metadata) == files_short[i]),]
    print(files_short[i])
    enoughdata = FALSE
    if (info$diagnosis == "Control" & info$protocol == "eyesclosed") {
      coll =  "black"
      count[1] = count[1] + 1
      if (count[1] > 2) enoughdata = TRUE
    } else if (info$diagnosis == "Epilepsy" & info$protocol == "eyesclosed") {
      coll = "grey"
      count[2] = count[2] + 1
      if (count[2] > 2) enoughdata = TRUE
    } else if (info$diagnosis == "Control" & info$protocol == "eyesopen") {
      coll = "darkblue"
      count[3] = count[3] + 1
      if (count[3] > 2) enoughdata = TRUE
    } else if (info$diagnosis == "Epilepsy" & info$protocol == "eyesopen") {
      coll = "blue"
      count[4] = count[4] + 1
      if (count[4] > 2) enoughdata = TRUE
    }
    if (enoughdata == FALSE) {
      x11()
      par(mfrow=c(7,2),oma=c(1,1,1,1),mar=c(1,1,1,1))
      psd = c()
      for (j in 1:14) {
        plot(myfun(x=data[1:t2,j]),type="l",ylim=yrange,col=coll)
        # plot(data[1:t2,j],type="l",ylim=yrange,col=coll)
        #           sp = spectrum(data[1:t2,j],plot=TRUE,col=coll)
        #           psd = rbind(psd,sp$spec)
      }
    }
    enoughdata = FALSE
  }
}


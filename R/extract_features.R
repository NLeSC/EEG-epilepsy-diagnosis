extract_features = function(cleandatadir,sf=128,n.levels=7,filtertypes,epochlength,fn){
  print("extract features")
  getfeatures = function(x,fns) {
    #--------------------------------------------------------
    # define generic function for feature extraction
    # Note: approx_entropy from pracma package is very slow
    # x: dataframe with multiple signals in columns over which the feature is calculated
    # fns: feature names that need to be calculated
    en.entropy = function(x) { # entropy
      en.entropy = entropy::entropy(entropy::discretize(as.numeric(x)+runif(n=length(x),min=0,max=1e-8)+0.0001,numBins=100))
    }
    pracma.samen = function(x) { #sample entropy
      pracma.samen = pracma::sample_entropy(x)
    }
    domfreq = function(x) {
      pp = stats::spectrum(x,plot=FALSE)
      domfreq = pp$freq[which.max(pp$spec)]
    }
    maxpow = function(x) {
      pp = stats::spectrum(x,plot=FALSE)
      maxpow = max(pp$spec^2)
    }
    RMS = function(x) {
      RMS = sqrt(mean(x^2))
    }
    zerocross = function(x) { # counts the number of zero crossing
      cc = rep(0,length(x))
      bf.fil = function(y) { # band-pass filter
        sf = 128
        lb = 0.2
        hb = 60
        Wc = matrix(0,2,1) 
        Wc[1,1] = lb / (sf/2) 
        Wc[2,1] = hb / (sf/2) 
        n = 4
        bf = signal::butter(n,Wc,type=c("pass")) 
        bf.fil = signal::filter(bf,y) 
      }
      x = bf.fil(x)
      cc[which(x >= 0)] = 1
      cc[which(x < 0)] = 0
      zerocross = length(which(abs(diff(cc)) == 1))
    }
    if (length(which(fns=="mean")) >0) A=aggregate(x=x,by=list(x$bands),mean)
    if (length(which(fns=="sd")) >0) A=aggregate(x=x,by=list(x$bands),sd)
    if (length(which(fns=="entropy")) >0) {
      A=aggregate(x=x,by=list(x$bands),en.entropy)
    }
    if (length(which(fns=="pracma.samen")) >0) {
      A=aggregate(x=x,by=list(x$bands),pracma.samen)
    }
    # if (length(which(fns =="RMS")) >0) A=aggregate(x=x,by=list(x$bands),RMS)
    return(A)
  }
  getid = function(x) {
    tmp = unlist(strsplit(x,"_"))[2]
    tmp2 = unlist(strsplit(tmp,"id"))[2]
    return(as.numeric(tmp2))
  }
  getprotocol = function(x) return(unlist(strsplit(x,"_"))[1])
  getdiagnosis = function(x) {
    tmp = unlist(strsplit(x,"_"))[5]
    tmp2 = unlist(strsplit(tmp,"gro"))[2]
    return(unlist(strsplit(tmp2,"[.]"))[1])
  }
  getdur = function(x) {
    tmp = unlist(strsplit(x,"_"))[3]
    tmp2 = unlist(strsplit(tmp,"dur"))[2]
    return(as.numeric(tmp2))
  }
  getepoch = function(x) {
    tmp = unlist(strsplit(x,"_"))[4]
    tmp2 = unlist(strsplit(tmp,"epoch"))[2]
    return(as.numeric(tmp2))
  }
  bf.fil = function(x,sf) { # low-pass filter
    x = x - mean(x)
#     hb = floor(sf/2) - 1
#     lb = 0.25 #we only have 4 seconds of data, so there is no point at looking at < 0.25 Hz
#     Wc = matrix(0,2,1)
#     Wc[1,1] = lb / (sf/2)
#     Wc[2,1] = hb / (sf/2)
#     bf = signal::butter(n=4,Wc,type=c("pass")) 
#     bf.fil = signal::filter(bf,x) 
      return(x)
  }
  #-------------------------------------------
  siglen = epochlength*sf
  # filenames
  files = list.files(cleandatadir,include.dirs=TRUE,full.names = TRUE)
  files_short = list.files(cleandatadir,include.dirs=FALSE,full.names = FALSE)
  # extract metadata from filenames
  metadata = data.frame(id = as.numeric(sapply(files_short,getid)),
                        protocol = as.character(sapply(files_short,getprotocol)),
                        diagnosis =  as.character(sapply(files_short,getdiagnosis)),
                        dur = as.numeric(sapply(files_short,getdur)),
                        epoch = as.numeric(sapply(files_short,getepoch)),row.names=files_short)
  metadata = metadata[order(metadata$id),]
  # initialize other parameters:
  S = c() #object in which features will be collected
  windowcnt = 0 #counter for total number of data windows (rows in output matrix)
  print(paste("filtertypes: ",paste(filtertypes,collapse=" ")))
  print(paste("N files: ",paste(length(files_short))))
  cnt = 0  #counter for showing process in console
  for (i in 1:length(files_short)) {#length(files_short)) { #unique id numbers
    cat(paste0(i," "))
    if (cnt == 10) { # print progress of processing
      prog = round((i/length(files_short)) * 1000)/ 10
      print(paste0(prog," %"))
      cnt = 0
    }
    cnt = cnt + 1
    cleaneddata = utils::read.csv(files[i]) # load files
    if (nrow(cleaneddata) >= siglen) { 
      Nwindows = floor(nrow(cleaneddata) / siglen)
      # cnt = 0
      for (windowi in 1:Nwindows) {
        windowcnt = windowcnt + 1
        windowdata = cleaneddata[(((windowi-1) * siglen) + 1):(siglen*windowi),2:15] # we are only interested in these columns
        windowdata = t(windowdata)
        wtdata = NULL
        windowdata = t(apply(windowdata,1,bf.fil,sf)) # subtract signal mean before performing wavelet analyses
        #------------------------------------------------------------
        # Wavelets:
        #1: 32-64 samplewindow #2: 16-32 samplewindow beta; #3: 8-16 samplewindow alpha;
        #4: 4-8 samplewindow theta #5: 2-4 samplewindow delta; #6: 1-2 samplewindow delta;    #7: 0.5-1 samplewindow delta
        # Sample frequency is 128Hz, so this roughly corresponds to:
        #1: 32-64 Hertz #2: 16-32 Hertz beta;    #3: 8-16 Hertz alpha;   #4: 4-8 Hertz theta
        #5: 2-4 Hertz delta; #6: 1-2 Hertz delta;    #7: 0.5-1 Hertz delta
        for (firi in filtertypes) { #filter types
          mymra = function(x){
            out = wavelets::mra(x,filter=firi, boundary="periodic",n.levels=n.levels)
            return(unlist(c(out@D)))
          }
          wtdata = t(apply(windowdata,1,mymra)) # apply multi-resolution analyses
          wtdata <- as.data.frame(wtdata) # these are the wavelets
          bands = c(rep(1:n.levels,each=siglen))  
          waveletdata = as.data.frame(t(wtdata))
          waveletdata$bands = bands
          #----------------------------------------
          # Connectivity on raw EEG data:
          pli_features = phaselagindex(EEGdata=t(windowdata),frequency=128)
          pli_features = as.data.frame(pli_features)
          names(pli_features) = paste0(names(pli_features),".raw.notapplicable") #,windowi)
          if (length(S) == 0) {
            S = pli_features
          } else {
            matchcolumns = which(names(S) %in% names(pli_features) == TRUE)
            if (length(matchcolumns) > 0) { #insert in matrix, rather than 
              S[windowcnt,matchcolumns] = pli_features
            } else {
              S = cbind(S,pli_features)
            }
          }
          # Connectivity on wavelets:
          t0 = Sys.time()
          for (ubi in unique(bands)) { 
            pli_features = phaselagindex(EEGdata=waveletdata[waveletdata$bands==ubi,which(names(waveletdata)!= "bands")],frequency=128)
            pli_features = as.data.frame(pli_features)
            names(pli_features) = paste0(names(pli_features),".",firi,".",ubi) #,".",windowcnti)
            matchcolumns = which(names(S) %in% names(pli_features) == TRUE)
            if (length(matchcolumns) > 0) { 
              S[windowcnt,matchcolumns] = pli_features
            } else {
              S = cbind(S,pli_features)
            }
          }
          # Other features (based on wavelets)
          for (featuresi in fn) { # features to summarize wavelets
            A = getfeatures(x=waveletdata,fns=featuresi)
            if (ncol(A) > 14) A = A[,2:15]
            tmp = as.numeric(apply(A,1,min)) # minimum features value across channels
            tmp = as.data.frame(t(tmp))
            tmp2 = as.numeric(apply(A,1,max)) # maximum features value across channels
            tmp2 = as.data.frame(t(tmp2))
            tmp3 = as.numeric(apply(A,1,mean)) # mean features value across channels
            tmp3 = as.data.frame(t(tmp3))
            tmp5 = as.numeric(apply(A,1,sd)) # sd features value across channelsspeed
            tmp5 = as.data.frame(t(tmp5))
            tmp4 = cbind(tmp,tmp2,tmp3,tmp5)
            colnames(tmp4) =c(paste0(1:n.levels,".",featuresi,".",firi,".min"),#,windowi),
                              paste0(1:n.levels,".",featuresi,".",firi,".max"), #,windowi),
                              paste0(1:n.levels,".",featuresi,".",firi,".mean"), #,windowi),
                              paste0(1:n.levels,".",featuresi,".",firi,".sd")) #,windowi)) # add labels and merge
            matchcolumns = which(names(S) %in% names(tmp4) == TRUE)
            if (length(matchcolumns) > 0) { #insert in matrix, rather than expand matrix size all the time
              S[windowcnt,matchcolumns] = tmp4
            } else {
              S = cbind(S,tmp4)
            }
            
          }
          rm(tmp4)
          S$window[windowcnt] = windowi
          S$fnames[windowcnt] = files_short[i]
        }
      }
    }
    
  }
  S = cbind(S[,which(names(S) != "window" & names(S) != "fnames")],S[,c("window","fnames")]) # update order of columns
  DAT = S
  LAB = metadata[which(rownames(metadata) %in% unique(S$fnames) == TRUE),]
  invisible(list(DAT=DAT,LAB=LAB))
  
}
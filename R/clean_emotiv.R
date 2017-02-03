clean_emotiv = function(datadir,metadatafile,outputdir,sf,gyrothreshold,
                        mindur,knownerrors,protocoltimes,referencegroup,condition_start_closed,
                        protocolvariable) {
  print("load and clean data")
  print(sf)
  #--------------------------------------------------------
  # subfunctions
  getfileinfo = function(datadir) {
    files = list.files(datadir,include.dirs=TRUE,full.names = TRUE)
    files_short = list.files(datadir,include.dirs=FALSE,full.names = FALSE)
    getid = function(x) {
      tmp = unlist(strsplit(x,"[.]csv"))[1]
      return(as.numeric(unlist(strsplit(tmp,"-"))[2]))
    }
    idnames = sapply(files_short,getid)
    
    fileinfo = data.frame(id=idnames,fnames_short=files_short,fnames_long=files,stringsAsFactors = FALSE)
    return(fileinfo)
  }
  definepoordata = function(eegdata,protocol) {
    # quality 0 = poor data
    # quality 1 =  healthy data
    # quality 2 = corrected artifacts
    poordataindices = which(eegdata$quality == 0 & eegdata$protocol == protocol) #only remove quality= 0 , but still use quality = 2 (corrected artifacts)
    startend = range(which(eegdata$protocol == protocol))
    poordataindices = c(startend[1],poordataindices,startend[2]) #add first and last sample
    return(poordataindices)
  }
  addqualityindicator= function(eegdata,knownerrors.df,gyrothreshold,id) {
    # add labels to unfit parts of the data based on qc scores and gyro:
    eegdata$quality = 1 # default is quality 1, which is good
    minqc = do.call(pmin,as.data.frame(eegdata[,22:36])) #minimum qc value per timestep across channels
    gyrox = abs(eegdata$GYROX-stats::median(eegdata$GYROX)) # absolute deviation from the median
    gyroy = abs(eegdata$GYROY-stats::median(eegdata$GYROY)) # absolute deviation from the median
    # Check for head movement: The unit of gyro is difficult to interpret
    # confirmed by http://www.bci2000.org/wiki/index.php/Contributions:Emotiv
    # no head movement seems to coincide with variations of around 4 units
    # so 30 units on a scale of >1000 would seem to at least ommit the extreme movement
    eegdata$quality[which(gyrox > gyrothreshold | gyroy > gyrothreshold | minqc < 3)] = 0 # data quality 0 is defined as poor
    #======================================================================
    ikne = which(knownerrors.df[,1] == id) # ikne = Id for the KNown Errors
    if (length(ikne) > 0) {
      for (j in 1:length(ikne)) { # loop over patient ids with known errors
        tmp1 = as.numeric(unlist(strsplit(knownerrors.df[ikne[j],2],":")))
        tmp2 = tmp1[1] * 60 + tmp1[2]
        tmp3 = as.numeric(unlist(strsplit(knownerrors.df[ikne[j],3],":")))
        tmp4 = tmp3[1] * 60 + tmp3[2]
        starti = tmp2*sf
        endi = tmp4*sf
        if (endi > nrow(eegdata)) endi = nrow(eegdata)
        if (starti > endi) {
          print(knownerrors.df[ikne[j],])
          warning(paste0("impossible start time "))
        }
        eegdata$quality[starti:endi] = 0 # data quality 0 is defined as poor
      }
    }
    return(eegdata)
  }
  addprotocollabels = function(protocoltimes,metadata,eegdata) {
    # adds protocol labels (eyes open or closed to the eegdata
    section1 = ((protocoltimes[1]+5)*sf):((protocoltimes[2]-5)*sf)
    endindex = min(c(nrow(eegdata),(protocoltimes[3]-5)*sf))
    section2 = (protocoltimes[2]*sf):endindex
    if (metadata[,protocolvariable][ind] == condition_start_closed) {
      closedi = section1; openi = section2
    } else {
      closedi = section2; openi = section1
    }
    eegdata$protocol = "unknown"
    eegdata$protocol[closedi] = "closed"
    eegdata$protocol[openi] = "open"
    return(eegdata)
  }
  write2file = function(x,outputdir,meta,dur,epoch,protocol) {
    # x = data to be saved
    utils::write.csv(x,paste0(outputdir,"/eyes",protocol,"_id",meta$subject.id,"_dur",
                              dur,"_epoch",epoch,"_gro",meta$Group,".csv"),row.names=FALSE)
    return(0)
  }
  correctartifact = function(x) {
    # identifies artificats and corrects them
    dx = diff(x)
    qt = quantile(abs(dx),probs=c(0.68),na.rm = TRUE) # assumption that at least 68% of data is not affected
    ww = which(abs(dx) > (5 * qt))
    dx[ww] = 0 #reset all differences larger than 5 sigma
    x = cumsum(c(x[1],dx))
    # apply rolling median function
    x = x - zoo::rollmedian(x,k=((sf*4)+1),align="center",
                            fill=c(stats::median(x[1:(sf*10)]),NA,stats::median(x[(length(x)-(sf*10)):length(x)])))
    return(x)
  }
  mymra = function(x){
    # apply multi resolution wavelet analyses
    out = wavelets::mra(x,filter="d6", boundary="periodic",n.levels=7)
    return(unlist(c(out@D)))
  }
  #--------------------------------------------------------
  # Start script
  extract_country = as.character(unlist(strsplit(outputdir,"/"))) # extract country from file directory specified by outputdir
  extract_country = extract_country[length(extract_country)]
  metadata = utils::read.csv(metadatafile) # get metadata
  fileinfo = getfileinfo(datadir) # extract id numbers from EEG filenames
  uid = sort(unique(fileinfo$id))
  metadata = merge(metadata,fileinfo,by.y="id",by.x="subject.id") # merge EEG fileinfo with metadata
  knownerrors.df = data.frame(matrix(unlist(knownerrors),ncol=3,byrow=T),stringsAsFactors = FALSE)
  amountdata = matrix(NA,length(uid),4) # initialize matrix to keep record of amount of data
  if (length(which(names(metadata) == "subject.id" | names(metadata) == "Group")) < 2) {
    warning('metadata needs to have subject.id and Group in variable names')
  }
  print(paste("N unique ids: ",paste(length(uid))))
  count_artificatcorrections = 0
  count_healthy = 0
  cnt = 0 #counter for showing analysis process in console
  correction_overview = matrix(NA,length(uid),14)
  for (i in uid) { # loop over unique id numbers derived from eeg files
    cat(paste0(i," "))
    ind2 = which(uid == i)  #index in the unique eegdata ids
    if (cnt == 10) {   # print progress after every 5 files
      prog = round((which(uid == i)/length(uid)) * 1000)/ 10
      print(paste0(prog," %"))
      cnt = 0
    }
    cnt = cnt + 1
    # get indices of matching information in metadata
    ind = which(metadata$subject.id == i)[1] #index metadata for this id
    if (length(ind) > 0) { # if there is metadata for this file
      amountdata[ind2,3] = i
      if (metadata$Group[ind] == referencegroup) {
        amountdata[ind2,4] = 0
      } else {
        amountdata[ind2,4] = 1
      }
      eegdata = read.csv(fileinfo$fnames_long[which(fileinfo$fnames_short == metadata$fnames_short[ind])])
      # ==============================================
      # Only include measurements with at least mindur minutes of data:
      if (nrow(eegdata) >= (mindur * 60 * sf)) {
        # add labels to unfit parts of the data based on researcher remarks:
        eegdata = addqualityindicator(eegdata,knownerrors.df,gyrothreshold,id=i)
        # add labels to measurement parts to clarify protocol, eyes open or closed:
        eegdata = addprotocollabels(protocoltimes,metadata,eegdata) 
        eegdata_corrected = apply(eegdata[,2:15],2,correctartifact)
        # # Old code: Option to keep track of when data was corrupted
        for (j in 1:14) {
          arti = which(abs(abs(eegdata[,j+1]-median(eegdata[,j+1])) - abs(eegdata_corrected[,j]))  > (4*sd(eegdata_corrected[,j])))
          if (length(arti) > 0) eegdata$quality[arti] = 2 # data quality 2 is considered an artifact
        }
        eegdata_raw = eegdata[,2:15]
        eegdata[,2:15] = eegdata_corrected #replace eeg data by corrected data
        # Now normalize the corrected eeg data
        eegdata[,2:15] = apply(eegdata[,2:15],2,FUN=function(x) x/sd(x))
        # also normalize the raw data
        eegdata_raw = apply(eegdata_raw,2,FUN=function(x) (x-median(x))/sd(x))
        # Wavelet extraction for the purpose of plotting, we do not use it here for feature extraction
        wtdata = t(apply(eegdata[,2:15],2,mymra)) # apply multi-resolution analyses
        waveletdata = as.data.frame(t(wtdata))
        siglen =nrow(eegdata)
        n.levels = 7
        bands = c(rep(1:n.levels,each=siglen))  
        waveletdata$bands = bands
        #=================================
        # Plotting
        createnewplot = TRUE # indicator of whether a new plot needs to be created
        YLIM = c(-6,16) # range for plot y-axis
        CX = 0.8
        lwdX = 0.3
        CL = c("red","#15353b","#005866","#46919d", "#ffda83","#ffb300")
        qc1line = rep(NA,nrow(eegdata))
        qc1line[which(eegdata$quality == 0)] = -2.5
        for (j in 1:14) {
          if (createnewplot == TRUE) {
            pngfile= as.character(unlist(strsplit(outputdir,"/")))
            pngfile = paste0(pngfile[1:length(pngfile)-1],collapse="/")
            pngfile = paste0(pngfile,"/images/",extract_country,"_rawdata_inspection_id",i,
                             "_file",ceiling(j/7),".pdf")
            pdf(file=pngfile,width = 7,height = 9)
            par(mfrow=c(7,2),mar=c(2,2,0,0),oma=rep(0,4),mgp=c(1.2,0.6,0))
            createnewplot = FALSE
          }
          if (j == 14) createnewplot = TRUE
          time  = (1:nrow(eegdata)) / sf
          if (j == 7 | j == 14) {
            XLAB = "time (sec)"
          } else {
            XLAB = ""
          }
          # first plot raw uncorrected eeg
          plot(time,eegdata_raw[,j],type="l",col=CL[1],ylim=YLIM,
               xlab=XLAB, ylab=expression(paste("normalized ",mu,"Volt")),
               cex.main=CX,lwd=lwdX,cex=CX,main="",
               cex.axis=CX,cex.lab=CX,bty="l",axes=FALSE,bty="l")
          if (j == 7 | j == 14) { # start new collumn
            axis(side = 1,at = seq(0,250,by=50),labels=seq(0,250,by=50),cex=CX,cex.axis=CX)
          }
          # axis(side = 2,at = seq(-5,5,by=500),labels=seq(-5,5,by=500),cex=CX,cex.axis=CX)
          # show corrrected EEG
          lines(time,eegdata[,j+1],type="l",lwd=lwdX,col="black")
          # show artifacts related to qc, movement or protocol compliance
          lines(time,qc1line,lwd=2,type="l",col=CL[4],lend=2)
          # show unexplained artifacts
          qc2line = rep(NA,nrow(eegdata))
          if (length(which(eegdata$quality == 2)) > 0) { # line to indicate regions with former artifacts
            
            qc2line[which(eegdata$quality == 2)] = -3.5 # highlight all parts 
            lines(time,qc2line,lwd=2,type="l",col=CL[6],lend=2) #
            correction_overview[ind2,j] = 1
          } else {
            correction_overview[ind2,j] = 0
          }
          # plot wavelets
          for (bandi in 1:n.levels) {
            wx = waveletdata[which(waveletdata$bands == bandi),j]
            nx = seq(0,250,by=50)
            lines(time,wx + (1*(bandi-1))+5,
                  type="l",col="darkgrey",lwd=lwdX,cex=CX,cex.axis=CX,cex.lab=CX)
          }
          if (j == 1) { # add legend only in the first of the 14 plots
            legend("topright",legend=c("raw EEG","corrected EEG","qc<3/gyro>30/compliance","artifacts","wavelets"),
                   col=c(CL[1],"black",CL[4],CL[6],"darkgrey"),
                   lwd=c(0.5,0.5,2,2,0.5),cex=0.5,ncol=3,lty=rep(1,5),bg="white",box.lwd=0.3)
          }
          text(x = 10,y = 5,labels = names(eegdata)[j+1],cex = 2)
          if (j == 14) dev.off() #| j ==7
        }
        #===================================
        # export longest continuous healthy protocol part to a csv
        for (protocol in c("open","closed")) {
          # get indices for poor data segments
          poordataindices = definepoordata(eegdata,protocol)
          # investigate how long are the healthy parts
          boutdur = diff(poordataindices) #lengths of good bouts
          # keep at least four seconds of data
          atleast4sec = which(boutdur > sf *4)
          if (length(atleast4sec) > 0) {
            boutdur_long = boutdur[atleast4sec]
            epoch = 1
            for (ii in 1:length(boutdur_long)) { # loop through all sufficiently long healthy bouts of data
              bi = which(boutdur == boutdur_long[ii])
              for (ci in bi) { #in case there are multiple bouts with the same length
                select = (poordataindices[ci]+1):(poordataindices[ci+1]-1)
                if (protocol == "open") {
                  amountdata[ind2,1] = epoch
                  dataopen = eegdata[select,]
                  write2file(x=dataopen,outputdir,meta=metadata[ind,],dur = floor(boutdur_long[ii]/sf),epoch,protocol) 
                } else if (protocol == "closed") {
                  amountdata[ind2,2] = epoch
                  dataclosed = eegdata[select,]
                  write2file(x=dataclosed,outputdir,meta=metadata[ind,],dur = floor(boutdur_long[ii]/sf),epoch,protocol) 
                }
                epoch = epoch + 1
              }
            }
          } else {
            print("not sufficient data") #no data is saved
          }
        }
        rm(eegdata)
      }
    }
  }
  invisible(list(amountdata=amountdata,correction_overview=correction_overview))
}
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
    poordataindices = which(eegdata$quality != 1 & eegdata$protocol == protocol)
    startend = range(which(eegdata$protocol == protocol))
    poordataindices = c(startend[1],poordataindices,startend[2]) #add first and last sample
    return(poordataindices)
  }
  addqualityindicator= function(eegdata,knownerrors.df,gyrothreshold) {
    eegdata$quality = 1 # eegdata$quality 1 is defined as good
    #======================================================================
    # add labels to unfit parts of the data based on qc scores and gyro:
    minqc = do.call(pmin,as.data.frame(eegdata[,22:36])) #minimum qc value per timestep across channels
    gyrox = abs(eegdata$GYROX-stats::median(eegdata$GYROX)) # absolute deviation from the median
    gyroy = abs(eegdata$GYROY-stats::median(eegdata$GYROY)) # absolute deviation from the median
    # Check for head movement: The unit of gyro is difficult to interpret
    # confirmed by http://www.bci2000.org/wiki/index.php/Contributions:Emotiv
    # no head movement seems to coincide with variations of around 4 units
    # so 30 units on a scale of >1000 would seem to at least ommit the extreme movement
    eegdata$quality[which(gyrox > gyrothreshold | gyroy > gyrothreshold | minqc < 3)] = 0 # data quality 0 is defined as poor
    #======================================================================
    
    ikne = which(knownerrors.df[,1] == i)
    if (length(ikne) > 0) {
      for (j in 1:length(ikne)) {
        tmp1 = as.numeric(unlist(strsplit(knownerrors.df[ikne[j],2],":")))
        tmp2 = tmp1[1] * 60 + tmp1[2]
        tmp3 = as.numeric(unlist(strsplit(knownerrors.df[ikne[j],3],":")))
        tmp4 = tmp3[1] * 60 + tmp3[2]
        if (((tmp4*sf)-(tmp2*sf)+1) > nrow(eegdata)) {
          print(tmp2)
          print(tmp4)
          print(nrow(eegdata))
        }
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
    section1 = ((protocoltimes[1]+5)*sf):((protocoltimes[2]-5)*sf)
    endindex = min(c(nrow(eegdata),(protocoltimes[3]-5)*sf))
    section2 = (protocoltimes[2]*sf):endindex
    # add labels
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
  #--------------------------------------------------------
  extract_country = as.character(unlist(strsplit(outputdir,"/")))
  extract_country = extract_country[length(extract_country)]
  metadata = utils::read.csv(metadatafile) # get metadata
  fileinfo = getfileinfo(datadir) # extract id numbers from filenames
  uid = sort(unique(fileinfo$id))
  metadata = merge(metadata,fileinfo,by.y="id",by.x="subject.id") # merge with metadata
  knownerrors.df = data.frame(matrix(unlist(knownerrors),ncol=3,byrow=T),stringsAsFactors = FALSE)
  cnt = 0 #counter for showing process in console
  amountdata = matrix(NA,length(uid),4) # initialize matrix to keep record of amount of data
  if (length(which(names(metadata) == "subject.id" | names(metadata) == "Group")) < 2) {
    warning('metadata needs to have subject.id and Group in variable names')
  }
  print(paste("N unique ids: ",paste(length(uid))))
  
  count_artificatcorrections = 0
  count_healthy = 0
  correction_overview = matrix(NA,length(uid),14)
  
  for (i in uid) { # loop over unique id numbers derived from eeg files
    
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
        eegdata = addqualityindicator(eegdata,knownerrors.df,gyrothreshold)
        # add labels to measurement parts to clarify protocol, eyes open or closed:
        eegdata = addprotocollabels(protocoltimes,metadata,eegdata) 
        #=======================================
        # Artifact correction
        removearti = function(x) {
          dx = diff(x)
          qt = quantile(abs(dx),probs=c(0.68),na.rm = TRUE) # assumption that at least 68% of data is not affected
          ww = which(abs(dx) > (5 * qt))
          dx[ww] = 0 #reset all differences larger than 25 sigma
          x = cumsum(c(x[1],dx))
          x = x - zoo::rollmedian(x,k=((sf*10)+1),align="center",
                                  fill=c(median(x[1:(sf*10)]),NA,median(x[(length(x)-(sf*10)):length(x)])))
          return(x)
        }
        eegdatabefore = eegdata
        windowdata = eegdata[,2:15]
        windowdata = apply(windowdata,2,removearti) # subtract mean
        mymra = function(x){
          out = wavelets::mra(x,filter="d6", boundary="periodic",n.levels=7)
          return(unlist(c(out@D)))
        }
        wtdata = t(apply(windowdata,2,mymra)) # apply multi-resolution analyses
        waveletdata = as.data.frame(t(wtdata))
        siglen =nrow(eegdata)
        n.levels = 7
        bands = c(rep(1:n.levels,each=siglen))  
        waveletdata$bands = bands
        wncount = 1
        for (j in 1:14) {
          eegdata[,j+1] = eegdata[,j+1] - median(eegdata[,j+1])
          YLIM = range(c(eegdata[,j+1]),na.rm=TRUE)
          YLIM[1] = YLIM[1] - (abs(YLIM[1]) * 0.3)
          YLIM[2] = YLIM[2] + (abs(YLIM[2]) * 0.3)
          if (wncount == 1) {
            pngfile= as.character(unlist(strsplit(outputdir,"/")))
            pngfile = paste0(pngfile[1:length(pngfile)-1],collapse="/")
            
            pngfile = paste0(pngfile,"/images/",extract_country,"_rawdata_inspection_id",i,"_file",ceiling(j/7),".png")
            png(file=pngfile,width = 7,height = 7,units="in",res=300)
            par(mfrow=c(7,2),mar=c(2,2,1.5,0),mgp=c(1,0.5,0),oma=rep(0,4))
            wncount = 2
          }
          if (j == 7) wncount = 1
          CX = 0.7
          lwdX = 0.4
          time  = (1:nrow(windowdata)) / sf
          plot(time,eegdata[,j+1],type="l",col="red",ylim=c(-600,850),xlab="time (sec)",ylab="microVolt",
               main=paste0("channel ",names(eegdata)[j+1]),cex.main=CX,lwd=lwdX,cex=CX,cex.axis=CX,cex.lab=CX,bty="l")
          lines(time,windowdata[,j],type="l",lwd=lwdX,col="black")
          qc1line = eegdata$quality
          qc1line[which(qc1line == 0)] = -525
          qc1line[which(qc1line == 1)] = NA # we do not want to highlight the good data
          lines(time,qc1line,lwd=2,type="l",col="blue",lend=2)
          qc2line = rep(NA,nrow(eegdata))
          arti = which(abs(abs(eegdata[,j+1]) - abs(windowdata[,j]))  > (2*sd(windowdata[,j]))) #| abs(eegdata[,j+1]) > 300
          if (length(arti) > 0) {
            qc2line[arti] = -525 # highlight all parts 
            lines(time,qc2line,lwd=2,type="l",col="green",lend=2)
            legend("topright",legend=c("raw EEG","corrected EEG","moving / qc<3 / protocol error","artifacts"),col=c("red","black","blue","green"),
                   lwd=c(0.5,0.5,2,2),cex=0.4,ncol=4,lty=rep(1,4),bg="white",box.lwd=0.3)
            eegdata$quality[arti] = 0
            correction_overview[ind2,j] = 1
          } else {
            correction_overview[ind2,j] = 0
          }
          for (bandi in 1:n.levels) {
            wx = waveletdata[which(waveletdata$bands == bandi),j]
            YLIM = c(-100,700) #range(wx)
            if (bandi == 1) {
            plot(time,wx,xlab="time (sec)",
                 type="l",col="black",
                 main=paste0("wavelets (short@bottom, long@top) from corrected EEG channel ",names(eegdata)[j+1]),
              ylim=YLIM,ylab="",axes=FALSE,#xlab="time",
                 cex.main=CX,lwd=lwdX,cex=CX,cex.axis=CX,cex.lab=CX)
              nx = seq(0,250,by=50)
              axis(1,at=nx,ps=nx,cex.axis=CX)
            } else {
              lines(time,wx + (100*(bandi-1)),xlab="time (sec)",
                   type="l",col="black",lwd=lwdX,cex=CX,cex.axis=CX,cex.lab=CX)
            }
          }
          if (j == 14 | j ==7) dev.off()
        }
        #===================================
        # investigate length of healthy time segments per protocol and
        # export longest continuous healthy protocol part to a csv
        for (protocol in c("open","closed")) {
          poordataindices = definepoordata(eegdata,protocol)
          boutdur = diff(poordataindices) #lengths of good bouts
          atleast4sec = which(boutdur > sf *mindur)
          if (length(atleast4sec) > 0) {
            # print("sufficient data") #no data is saved
            boutdur_long = boutdur[atleast4sec]
            epoch = 1
            for (ii in 1:length(boutdur_long)) { # loop through all long bouts of healthy data
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
                  # remove artifacts and keep track of when this is done
                  write2file(x=dataclosed,outputdir,meta=metadata[ind,],dur = floor(boutdur_long[ii]/sf),epoch,protocol) 
                  
                }
                epoch = epoch + 1
              }
            }
          } else { #store if there is at least 4 seconds of data
            print("not sufficient data") #no data is saved
          }
        }
        rm(eegdata)
      }
    }
  }
  invisible(list(amountdata=amountdata,correction_overview=correction_overview))
}
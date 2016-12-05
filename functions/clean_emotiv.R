clean_emotiv = function(datadir,metadatafile,outputdir,sf,gyrothreshold,mindur,knownerrors,protocoltimes,referencegroup,condition_start_closed) {
  print("load and clean data")
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
    gyrox = abs(eegdata$GYROX-median(eegdata$GYROX)) # absolute deviation from the median
    gyroy = abs(eegdata$GYROY-median(eegdata$GYROY)) # absolute deviation from the median
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
#         print(paste0(length(eegdata$quality)," ",eegdata$quality[1]," ",knownerrors.df[ikne[j],1]," ",
#                      knownerrors.df[ikne[j],2]," ",knownerrors.df[ikne[j],3]))
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
    write.csv(x,paste0(outputdir,"/eyes",protocol,"_id",meta$subject.id,"_dur",
                       dur,"_epoch",epoch,"_gro",meta$Group,".csv"),row.names=FALSE)
    return(0)
  }
  #--------------------------------------------------------
  metadata = read.csv(metadatafile) # get metadata
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
      # Only include measurements with at least mindur minutes of data:
      if (nrow(eegdata) >= (mindur * 60 * sf)) {
        # add labels to unfit parts of the data based on researcher remarks:
        eegdata = addqualityindicator(eegdata,knownerrors.df,gyrothreshold)
        # add labels to measurement parts to clarify protocol, eyes open or closed:
        eegdata = addprotocollabels(protocoltimes,metadata,eegdata) 
        # investigate length of healthy time segments per protocol and
        # export longest continuous healthy protocol part to a csv
        for (protocol in c("open","closed")) {
          poordataindices = definepoordata(eegdata,protocol)
          boutdur = diff(poordataindices) #lengths of good bouts
          atleast4sec = which(boutdur > sf *mindur)
          if (length(atleast4sec) > 0) {
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
  return(amountdata)
}
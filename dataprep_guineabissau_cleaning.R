rm(list=ls())
graphics.off()
setwd("/home/vincent/utrecht/EEG-epilepsy-diagnosis")
datadir =  "/media/windows-share/EEGs_Guinea-Bissau__16-06-2016" #"data/eeg"
metadatafile = "/media/windows-share/subject.id_with_meta-info__anonymized.csv"
outputdir =  "/media/windows-share/EEGs_Guinea-Bissau_cleaned" 
sf = 128 #sample frequency
gyrothreshold = 30 #gyro unit deviations from the median that are not tolerated
#--------------------------------------------------------
# get metadata
metadata = read.csv(metadatafile) # read.csv("data/labels.csv")
# extract id numbers from filenames
files = list.files(datadir,include.dirs=TRUE,full.names = TRUE)
files_short = list.files(datadir,include.dirs=FALSE,full.names = FALSE)
getid = function(x) {
  tmp = unlist(strsplit(x,"[.]csv"))[1]
  return(as.numeric(unlist(strsplit(tmp,"-"))[2]))
}
idnames = sapply(files_short,getid)
uid = sort(unique(idnames))
fileinfo = data.frame(id=idnames,fnames=files_short)
# merge with metadata
metadata = merge(metadata,fileinfo,by.y="id",by.x="subject.id") 
# define known errors based on Research Remarks (this is hardcoded, but could be extracted from a file in the future)
knownerrors = list(c(10,"4:30","5:00"), # => id, starttime, endtime of the problematic period
                   c(13,"3:00","3:30"), # if there are more periods per id then use new entry
                   c(20,"0:55","1:10"),
                   c(22,"2:55","3:10"),
                   c(29,"3:00","3:30"),
                   c(40,"3:00","3:30"),
                   c(44,"3:00","4:00"),
                   c(58,"3:00","3:30"),
                   c(61,"3:00","5:00"),
                   c(67,"3:00","3:30"),
                   c(84,"2:20","2:30"),
                   c(92,"2:10","2:20"))
knownerrors.df = data.frame(matrix(unlist(knownerrors),ncol=3,byrow=T),stringsAsFactors = FALSE)
# initialize some parameters before loading and processing all the data
cnt = 0 #counter for showing process in console
print(paste("N unique ids: ",paste(length(uid))))
amountdata = matrix(NA,length(uid),3)
for (i in uid) { #loop over unique id numbers
  
  if (cnt == 5) { # print progress after every 5 files
    prog = round((i/length(uid)) * 1000)/ 10
    print(paste0(prog," %"))
    cnt = 0
  }
  cnt = cnt + 1
  ind = which(metadata$subject.id == i)[1] #which metadata belongs to this id
  if (length(ind) > 0) { # if there is metadata for this file
    amountdata[ind,3] = i
    D = read.csv(files[which(files_short == metadata$fnames[ind])])
    # we do expect measurements of 5 minutes, but lets include all with > 4 minutes for now:
    if (nrow(D) >= (4 * 60 * sf)) {
      D$quality = 1 # D$quality 1 is defined as good
      #======================================================================
      # add labels to unfit parts of the data based on qc scores and gyro:
      minqc = do.call(pmin,as.data.frame(D[,22:36])) #minimum qc value per timestep across channels
      gyrox = abs(D$GYROX-median(D$GYROX)) # absolute deviation from the median
      gyroy = abs(D$GYROY-median(D$GYROY)) # absolute deviation from the median
      # Check for head movement: The unit of gyro is difficult to interpret
      # confirmed by http://www.bci2000.org/wiki/index.php/Contributions:Emotiv
      # no head movement seems to coincide with variations of around 4 units
      # so 30 units on a scale of >1000 would seem to at least ommit the extreme movement
      D$quality[which(gyrox > gyrothreshold | gyroy > gyrothreshold | minqc < 3)] = 0 # data quality 0 is defined as poor
      #======================================================================
      # add labels to unfit parts of the data based on researcher remarks:
      ikne = which(knownerrors.df[,1] == i)
      if (length(ikne) > 0) {
        for (j in 1:length(ikne)) {
          tmp1 = as.numeric(unlist(strsplit(knownerrors.df[ikne[j],2],":")))
          tmp2 = tmp1[1] * 60 + tmp1[2]
          tmp3 = as.numeric(unlist(strsplit(knownerrors.df[ikne[j],3],":")))
          tmp4 = tmp3[1] * 60 + tmp3[2]
          D$quality[(tmp2*sf):(tmp4*sf)] = 0 # data quality 0 is defined as poor
        }
      }
      #======================================================================
      # Add labels to measurement parts to clarify protocol, eyes open or closed:
      # identify the two measurement sections
      section1 = (125*sf):(235*sf)
      endindex = min(c(nrow(D),355*sf))
      section2 = (245*sf):endindex #nrow(D)
      # add labels
      if (metadata$Eyes.condition[ind] == "closed-3min-then-open-2min") {
        closedi = section1; openi = section2
      } else {
        closedi = section2; openi = section1
      
      }
      D$protocol = "unknown"
      D$protocol[closedi] = "closed"
      D$protocol[openi] = "open"
      #======================================================================
      # investigate length of healthy time segments per protocol and
      # export longest continuous healthy protocol part to a csv
      for (protocol in c("open","closed")) {
        poordataindices = which(D$quality != 1 & D$protocol == protocol)
        startend = range(which(D$protocol == protocol))
        poordataindices = c(startend[1],poordataindices,startend[2]) #add first and last sample
        boutdur = diff(poordataindices) #lengths of good bouts
        atleast4sec = which(boutdur > sf *4)
        if (length(atleast4sec) > 0) {
          boutdur_long = boutdur[atleast4sec]
          print(paste0("N blocks: ",length(boutdur_long)))
          for (ii in 1:length(boutdur_long)) {
            bi = which(boutdur == boutdur_long[ii])
            select = (poordataindices[bi]+1):(poordataindices[bi+1]-1)
            if (protocol == "open") {
              dataopen = D[select,]
              amountdata[ind,1] = boutdur_long[ii]/sf
              write.csv(dataopen,paste0(outputdir,"/eyesopen_id",metadata$subject.id[ind],"_dur",
                                        floor(boutdur_long[ii]/sf),
                                        "_epoch",ii,"_gro",metadata$Group[ind],".csv"),row.names=FALSE)
            }
            if (protocol == "closed") {
              dataclosed = D[select,]
              amountdata[ind,2] = boutdur_long[ii]/sf
              write.csv(dataclosed,paste0(outputdir,"/eyesclosed_id",metadata$subject.id[ind],"_dur",
                                          floor(boutdur_long[ii]/sf),"_gro",metadata$Group[ind],".csv"),row.names=FALSE)
            }
          }
        } else { #store if there is at least 4 seconds of data
          print("not sufficient data") #no data is saved
        }
      }
      rm(D)
    }
  }
}
print(paste0("successful open: ",length(which(is.na(amountdata[,1]) == FALSE)) / nrow(amountdata)))
print(paste0("successful closed: ",length(which(is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
print(paste0("succesful both: ",length(which(is.na(amountdata[,1]) == FALSE &
                                               is.na(amountdata[,2]) == FALSE)) / nrow(amountdata)))
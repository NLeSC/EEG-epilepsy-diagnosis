# script to convert all output into tables and figures for the paper
rm(list=ls())
graphics.off()
path = "/media/windows-share/EEG/"



# Figure 1:
correctartifact = function(x) {
  sf = 128
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
datadir = "/media/windows-share/EEG/EEGs_Nigeria"
file2plot = dir(datadir,full.names = TRUE)
eegdata = read.csv(file2plot[5])

eegdata = eegdata[(30*128):(90*128),]
S = eegdata[1:((nrow(eegdata)/128)*128),]
time = seq(1/128,nrow(S)/128,by=1/128)

SAF4_corrected = correctartifact(S$AF4)
outfile = paste0(path,paste0("/Figure1.jpeg"))
jpeg(filename=outfile, units="in",width = 7,height= 3.5,res=600,pointsize = 12)
par(mar=c(4,4.5,3,1),mfrow=c(1,2))
YLIM = c(-400,400)
plot(time,S$AF4,type="l",xlab="time (seconds)",ylab=expression(paste(mu,"Volt")),
     cex=0.2,bty="l",cex.axis=0.8,cex.lab=1,main="original",ylim=YLIM)
plot(time,SAF4_corrected,type="l",xlab="time (seconds)",ylab=expression(paste(mu,"Volt")),
     cex=0.2,bty="l",cex.axis=0.8,cex.lab=1,main="corrected",ylim=YLIM)
# plot(time,S$F4,type="l",xlab="",ylab=expression(paste(mu,"Volt F4")),cex=0.5,bty="l")
# plot(time,S$FC5,type="l",xlab="time (seconds)",ylab=expression(paste(mu,"Volt FC5")),cex=0.5,bty="l")
dev.off()


# eegdata_corrected = apply(eegdata[,2:15],2,correctartifact)
#======================================
# count how much data is available
winsize = 4
# for (country in 1:2) {
#   if (country == 1) {
#     cat("Nigeria:\n")
#     path_cleaned = "/media/windows-share/EEG/EEGs_ni_cleaned"
#     contr = "control"; epi = "epilepsy"
#     
#   } else {
#     cat("Guinea-Bissau:\n")
#     path_cleaned = "/media/windows-share/EEG/EEGs_gb_cleaned"
#     contr = "Control"; epi = "Epilepsy"
#   }
#   fnames = dir(path_cleaned)
#   getinfo = function(x) {
#     tmp = as.character(unlist(strsplit(x,"_")))
#     protocol = tmp[1]
#     id = unlist(strsplit(tmp[2],"d"))[2]
#     dur= as.numeric(unlist(strsplit(tmp[3],"r"))[2])
#     epoch= as.numeric(unlist(strsplit(tmp[4],"h"))[2])
#     tmp2 = unlist(strsplit(tmp[5],"gro"))[2]
#     diagnosis= as.character(unlist(strsplit(tmp2,"[.]"))[1])
#     return(list(protocol=protocol,id=id,dur=dur,epoch=epoch,diagnosis=diagnosis))
#   }
#   ad = as.data.frame(t(sapply(fnames,getinfo)))
#   ad$protocol = as.factor(as.character(ad$protocol))
#   ad$diagnosis = as.factor(as.character(ad$diagnosis))
#   ad$id = as.numeric(ad$id)
#   ad$dur = as.numeric(ad$dur)
#   ad$epoch = as.numeric(ad$epoch)
#   ad$dur2 = floor(as.numeric(ad$dur) / winsize)
#   ad = aggregate(. ~ protocol + id + diagnosis,data=ad,FUN=sum)
#   for (prt in c("eyesclosed","eyesopen")) {
#     for (diag in c(contr,epi)) {
#       select = which(ad$protocol == prt & ad$diagnosis == diag & ad$dur2 >= 1)
#       if (length(select) > 0) {
#         cat(paste0(prt," ",diag," ",length(unique(ad[select,]$id))," with on averageld ",
#                    round(mean(ad[select,]$dur2),digits=1)," windows per persoon\n"))
#         
#       } else {
#         print(paste0(prt," ",diag," ",0))
#       }
#     }
#   }
# }
# 
# # Count how often correction has been applied:
# print(paste0("epoch=",winsize," value=",0))
# load(paste0(path ,"EEgs_logs/correctionoverview_gb.RData"))
# print(paste0("GB rowmeans=0: ",length(which(rowMeans(correction_overview) == 0))))
# print(paste0("GB rowsum<3: ",length(which(rowSums(correction_overview) < 3))))
# print(paste0("GB files: ",length(which(is.na(rowMeans(correction_overview)) == FALSE))))
# load(paste0(path,"EEgs_logs/correctionoverview_ni.RData"))
# print(paste0("NI rowmeans=0: ",length(which(rowMeans(correction_overview) == 0))))
# print(paste0("NI rowsum<3: ",length(which(rowSums(correction_overview) < 3))))
# print(paste0("NI files: ",length(which(is.na(rowMeans(correction_overview)) == FALSE))))

# kk
# TABLE 2
limit2sdfeatutes = TRUE
for (winsize in c(4)) {
  NCO = 8 # number of conditions
  output = data.frame(protocol=rep(0,NCO ),countrytrain=rep(" ",NCO ),evalcountry=rep("d",NCO ),
                      aggperid=rep("d",NCO),winningwavelet=rep("d",NCO),
                      sens=rep("d",NCO ),kappa=rep("d",NCO ),auc=rep("d",NCO ),acc=rep("d",NCO ),stringsAsFactors = FALSE)
  cnt = 1  
  for (countrytrain in c("gb","ni")) {
    for (evalcountry in c("gb","ni")) {
      # evalcountry = countrytrain
      for (aggperid in  c(TRUE,FALSE)) { #,FALSE
        for (proto_i in 2) {
          file = paste0("EEGs_evaluation/seedcomparison_",proto_i,"_dur",winsize,"_country",countrytrain,"_perid",
                        aggperid,"_evalcountry",evalcountry,"_SDonly",limit2sdfeatutes,".RData")
          load(paste0(path,file))
          M = lapply(evse,mean)
          S = lapply(evse,sd)
          output$protocol[cnt] = proto_i
          output$countrytrain[cnt] = countrytrain
          output$evalcountry[cnt] = evalcountry
          output$aggperid[cnt] = aggperid
          output$winningwavelet[cnt] = names(which.max(table(evse$winningmodel)))[1] #as.character(evse$winningmodel[round(nrow(evse)/2)])
          output$sens[cnt] = paste0(round(M$test.sens*100,digits=0),"se",round((S$test.sens/sqrt(19))*100,digits=0))
          output$acc[cnt] = paste0(round(M$test.acc*100,digits=0),"se",round((S$test.acc/sqrt(19))*100,digits=0))
          output$kappa[cnt] = paste0(round(M$test.kappa,digits=2),"se",round(S$test.kappa/sqrt(19),digits=2))
          output$auc[cnt] = paste0(round(M$test.auc,digits=2),"se",round(S$test.auc/sqrt(19),digits=2))
          # print(as.character(evse$winningmodel[round(nrow(evse)/2)]))
          cnt = cnt + 1
        }
      }
    }
  }
  write.csv(output,file=paste0(path,"table2_",winsize,".csv"))
  print(output)
}


library(randomForest)
if (limit2sdfeatutes == TRUE) {
  Nfeatures = 28
} else {
  Nfeatures = 140
}
for (winsize in c(4)) { #,4
  for (aggperid in c(FALSE,TRUE)) {
    outfile = paste0(path,paste0("/variablecontributions_aggregation",aggperid,".jpeg"))
    if (limit2sdfeatutes == TRUE) {
      jpeg(filename=outfile, units="in",width = 7,height= 5,res=600,pointsize = 10)
      par(oma=c(0,0,0,0),mar=c(5.5,5,3,1),mfrow=c(1,2))
    } else {
      jpeg(filename=outfile, units="in",width = 7,height= 7,res=600,pointsize = 10)
      par(oma=c(0,0,0,0),mar=c(1,1,0,1),mfrow=c(1,2))
    }
    for (countrytrain in c("gb","ni")) {
      for (proto_i in 2) {
        # outfile = paste0(path,paste0("/variablecontributions_",aggperid,"_",countrytrain," ",proto_i,".jpeg"))
        
        
        VI = data.frame(x1=rep(0,Nfeatures),x2=rep(0,Nfeatures),x3=rep(0,Nfeatures),x4=rep(0,Nfeatures))
        cnt = 1
        filebeginning = paste0("bestmodel_",proto_i,"_dur",winsize,"_country",countrytrain,
                               "_perid",aggperid,"_SDonly",limit2sdfeatutes)
        fnames = dir(paste0(path,"EEGs_bestmodels"))
        fnames2 = fnames[which(unlist(lapply(fnames,function(x) {length(unlist(strsplit(x,"_seed")))})) == 2)]
        fnames3 = unlist(lapply(fnames2,function(x) {unlist(strsplit(x,"_seed"))[1]}))
        fnames4 = fnames[which(fnames3 %in% filebeginning == TRUE)]
        # kk
        for (jj in 1:length(fnames4)) {
          file = paste0(path,"EEGs_bestmodels/",fnames4[jj])
          load(file)
          varimp = importance(best_model$finalModel,type=2,scale=FALSE)
          VI[,cnt] = varimp / sum(varimp)
          rownames(VI) = rownames(varimp)
          cnt = cnt +1
        }
        VI_mean = rowMeans(VI)
        rnames = rownames(VI)
        
        # extract wavelet level from variable names non pli, for non-wavelet variables assigned level 0
        level = as.numeric(sapply(rnames,function(x) {unlist(strsplit(x,"[.]"))[1]}))
        level[which(is.na(level) == TRUE)] = 0
        # construct wavelet variable names without wavelet level
        wvarnames_without_level = sapply(rnames[which(level != 0)],function(x) {
          b = unlist(strsplit(x,"[.]")); return(paste0(b[2:length(b)],collapse="."))
        })
        # extract wavelet level from variable names from pli, for non-wavelet variables assigned level 0
        plilevel = as.character(sapply(rnames[which(level == 0)],function(x) {
          x = unlist(strsplit(x,"[.]")); return(x[length(x)])
        }))
        pli = rep(0,length(level))
        pli[which(level == 0)] = 1
        level[which(level == 0)[which(plilevel != "notapplicable")]] = as.numeric(plilevel[which(plilevel != "notapplicable")])
        
        wvarnames_with_level = as.character(sapply(rnames[which(level != 0 & pli == 1)],function(x) {
          b = unlist(strsplit(x,"[.]")); return(paste0(b[1:(length(b)-1)],collapse="."))
        }))
        L0 = which(level == 0)
        if (limit2sdfeatutes == FALSE) {
          NV = length(unique(wvarnames_without_level)) + length(L0) + 14
        } else {
          NV = length(unique(wvarnames_without_level))
        }
        

        VI_matrix = matrix(0,NV,7)
        if (limit2sdfeatutes == FALSE) {
          VI_matrix[1:length(L0),1] = VI_mean[L0] # raw pli
        }
        for (i in 1:7) {
          VI_matrix[(length(L0)+1):NV,i] = VI_mean[which(level == i)]
        }
        
        dfrownames = c(rnames[L0],unique(wvarnames_with_level),unique(wvarnames_without_level))
        rm.notapp = function(x) {
          if (length(unlist(strsplit(x,".notapp"))) > 1) {
            x = unlist(strsplit(x,".notapp"))[1]
          }
          return(x)
        }
        replace.d = function(x) {
          xt = paste0("ppp",x,"ppp")
          for (ww in c(".d10",".d6",".d2")) {
            if (length(unlist(strsplit(xt,ww))) > 1) {
              xt = unlist(strsplit(xt,ww))
              xt = paste0(xt[1],".wavelet",xt[2])
              x = paste0(unlist(strsplit(xt,"ppp")),collapse = "")
            }
          }
          return(x)
        }
        dfrownames = as.character(sapply(dfrownames,replace.d))
        dfrownames = sapply(dfrownames,rm.notapp)
        dfrownames = gsub("ninetyp", "p90", dfrownames)
        dfrownames = gsub("pli", "_pli", dfrownames)
        reorder.sd.wavelet = function(x) {
          xt = paste0("ppp",x,"ppp")
          if (length(unlist(strsplit(xt,"d.wavelet."))) > 1) {
            xt = unlist(strsplit(xt,"d.wavelet."))
            xt = paste0(xt[2],".sd.wavelet") #,xt[2])
            x = paste0(unlist(strsplit(xt,"ppp")),collapse = "")
          }
          return(x)
        }
        dfrownames = sapply(dfrownames,reorder.sd.wavelet )
        df = data.frame(VI_matrix,row.names=dfrownames)
        names(df) = paste0("Level",1:7)
        # library(corrplot)

        TTL_country = "Guinea Bissau"
        if (countrytrain == "ni") TTL_country = "Nigeria"

        plot(as.numeric(df[1,]),type="p",pch=0,main=TTL_country,ylim=c(0,0.15),cex.lab=1.1,
             bty="l",axes = FALSE,xlab="Wavelet",ylab="Importance (mean decrease in Gini)",mgp=c(3.5,1,0)) #,cex.axis=0.4)
        axis(1,labels=c("Gamma","Beta","Alpha","Theta","Delta 2-4Hz","Delta 1-2Hz","Delta <1Hz"),
             at = 1:7,las=3,cex.axis=0.8)
        xx = seq(0,0.15,by=0.05)
        axis(2,labels=paste0(xx),at = xx,cex.axis=0.8)
        lines(as.numeric(df[2,]),type="p",pch=1)
        lines(as.numeric(df[3,]),type="p",pch=2)
        lines(as.numeric(df[4,]),type="p",pch=3)
        legend("topleft",legend = c("minimum sd","maximum sd","mean sd","sd of the sd"),pch = c(0:3,19),bty="o",cex=0.9) #rownames(df)
        # TTL = paste0(TTL_country)
        # corrplot(t(df),title = TTL,is.corr=FALSE,na.label = " ", cl.lim=c(0,0.18),cl.ratio=0.02,
        #          method="color",number.cex=0.9,cl.cex=1,tl.cex=0.9,tl.col="black",mar=c(0,0,2,0)) #circle"
        
        
      }
    }
    dev.off()
  }
}

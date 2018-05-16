# script to convert all output into tables and figures for the paper
# This script is not intended to be a generic function, but more a log of the code we used to generate some
# of the (prelimenary) tables and figures
rm(list=ls())
graphics.off()
path = "/media/vincent/EEG/"


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
datadir = "/media/vincent/EEG/EEGs_Nigeria"
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
  
  # outfile = paste0(path,paste0("/variablecontributions_aggregation",aggperid,".jpeg"))
  outfile = paste0(path,paste0("/variablecontributions.jpeg"))
  if (limit2sdfeatutes == TRUE) {
    jpeg(filename=outfile, units="in",width = 7,height= 7,res=600,pointsize = 10)
    par(oma=c(0,0,0,0),mar=c(3.5,3,3,0),mfrow=c(2,2))
  } else {
    jpeg(filename=outfile, units="in",width = 7,height= 7,res=600,pointsize = 10)
    par(oma=c(0,0,0,0),mar=c(1,1,0,1),mfrow=c(1,2))
  }
  for (aggperid in c(FALSE,TRUE)) {
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
        confinter = function(x) {
          x = (sd(x) / sqrt(length(x))) * 1.96
        }
        
        VI_sd = apply(VI,1,confinter)
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
        VIsd_matrix = matrix(0,NV,7)
        if (limit2sdfeatutes == FALSE) {
          VI_matrix[1:length(L0),1] = VI_mean[L0] # raw pli
          VIsd_matrix[1:length(L0),1] = VI_sd[L0] # raw pli
        }
        for (i in 1:7) {
          VI_matrix[(length(L0)+1):NV,i] = VI_mean[which(level == i)]
          VIsd_matrix[(length(L0)+1):NV,i] = VI_sd[which(level == i)]
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
        
        dfsd = data.frame(VIsd_matrix,row.names=dfrownames)
        names(dfsd) = paste0("Level",1:7)
        
        # library(corrplot)
        if (aggperid == TRUE) {
          maxy = 0.22
          heigthlegend = 0.2
        } else {
          maxy = 0.22
          heigthlegend = 0.2 #0.015
        }
        xspacing = 0.15
        CX = 1.1
        mgpVAL = c(2,0.5,0)
        
        TTL_country = "Guinea Bissau"
        if (countrytrain == "ni") TTL_country = "Nigeria"
        
        if (aggperid == TRUE) {
          TTL_country = paste0(TTL_country," (aggregated)")
        } else {
          TTL_country = paste0(TTL_country," (not aggregated)")
        }
        
        
        xposition = (1:7) -0.5 + (xspacing*(-1.5))
        if (countrytrain != "ni") {
          plot(xposition,as.numeric(df[1,]),type="p",pch=15,main=TTL_country,ylim=c(0,maxy),xlim=c(0,(7+(0.15*(4-3.5)))),cex.lab=1,cex.main=1.2,cex=CX,
               bty="n",axes = FALSE,xlab="Wavelet",ylab="Importance (decrease in Gini)",mgp=mgpVAL) #,cex.axis=0.8)
        } else {
          plot(xposition,as.numeric(df[1,]),type="p",pch=15,main=TTL_country,ylim=c(0,maxy),xlim=c(0,(7+(0.15*(4-3.5)))),cex.lab=1,cex.main=1.2,cex=CX,
               bty="n",axes = FALSE,xlab="Wavelet",ylab="",mgp=mgpVAL) #,cex.axis=0.8)
        }
        
        axis(1,labels=c("Gamma","Beta","Alpha","Theta","Delta1","Delta2","Delta3"),
             at = (1:7)-0.5,las=1,cex.axis=0.6)
        xx = seq(0,maxy,by=0.04)
        if (countrytrain != "ni") axis(2,labels=paste0(xx),at = xx,cex.axis=0.7)
        xposition = (1:7) -0.5+ (xspacing*(-0.5))
        lines(xposition,as.numeric(df[2,]),type="p",pch=16,cex=CX)
        xposition = (1:7) -0.5+ (xspacing*(0.5)) 
        lines(xposition,as.numeric(df[3,]),type="p",pch=17,cex=CX)
        xposition = (1:7) -0.5+ (xspacing*(1.5)) 
        lines(xposition,as.numeric(df[4,]),type="p",pch=18,cex=CX+0.2)
        for (ai in 1:4) {
          xposition = (1:7) -0.5+ (xspacing*(ai-2.5)) 
          arrows(xposition, as.numeric(df[ai,]-dfsd[ai,]), xposition, as.numeric(df[ai,]+dfsd[ai,]), length=0.02, angle=90, code=3)
        }
        
        if (countrytrain == "ni") legend(x = 0.8,y = heigthlegend,legend = c("minimum","maximum","mean","standard deviation"),pch = c(15:18),bty="o",cex=0.9,pt.cex = c(0.8,0.8,0.8,1),ncol=2) #rownames(df)
        figure3data = rbind(df,dfsd)
        write.csv(figure3data,file=paste0(path,"/figure3data_aggregated",aggperid,".csv"))
        # TTL = paste0(TTL_country)
        # corrplot(t(df),title = TTL,is.corr=FALSE,na.label = " ", cl.lim=c(0,0.18),cl.ratio=0.02,
        #          method="color",number.cex=0.9,cl.cex=1,tl.cex=0.9,tl.col="black",mar=c(0,0,2,0)) #circle"
        
        
      }
    }
    
  }
  dev.off()
}

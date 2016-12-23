rm(list=ls())
graphics.off()
library(emotivepilepsy)

setwd("/home/vincent/utrecht")
shareddrive = "/media/windows-share/EEG"

funcfiles = list.files("emotivepilepsy/R",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) source(i)

trainbestmodel = FALSE

for (aggregateperid in c(TRUE)) {
  for (proto_i in  2:1) { #"open" =1 #closed= 2
    # proto_i = 1  # "open" 1 #closed 2
    logfile = "data/log_guinneabissau.csv" # not used when uselog = FALSE
    load(file=paste0(shareddrive,"/features_and_bestmodels/features_nigeria_10.RData")); logdur = 10
    LAB$diagnosis = as.character(LAB$diagnosis)
    LAB$diagnosis[which(LAB$diagnosis == "control")] = "Control"
    LAB$diagnosis[which(LAB$diagnosis == "epilepsy")] = "Epilepsy"
    LAB$diagnosis = as.factor(LAB$diagnosis)
    
    # aggregateperid = FALSE
    print("============================")
    print(proto_i)
    print(aggregateperid)
    
    
    # tidy up formatting to be suitable for classifier training
    RDL = reformat_DATLAB(DAT,LAB,aggregateperid=aggregateperid) # aggregate per unique id
    DAT =RDL$DAT
    LAB = RDL$LAB
    
    #===============================================================
    # split data in training, validation and test set
    P = split_data(LAB,DAT,logfile = logfile,proto_i=proto_i,split=c(20,20),
                   uselog = FALSE,logdur=logdur)
    LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
    # generate dictionary of model characteristics
    modeldict = create_modeldict(DAT)
    #===============================================================
    # train models or loaded previously trained model
    if (trainbestmodel == TRUE) {
      trainingresults = train_model(DATtrain,LABtrain,DATval,LABval,modeldict,classifier = "rf") #,aggregateperid=aggregateperid
      best_model = trainingresults$best_model
      modelcomparison = trainingresults$result
      fes = trainingresults$fes
      country = "ni"
      bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",aggregateperid,".RData")
      save(best_model,fes,file=bestmodelfile) # Save best model
      rm(best_model,fes)
    } else {
      country = "gb" # load model from other country
      bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",aggregateperid,".RData")
      # LABtest = rbind(LABval) # ignore test data, and only evaluate on training and validation data
      # DATtest = rbind(DATval)
      LABtest = LABval # LABtest = rbind(LABval) # ignore test data, and only evaluate on training and validation data
      DATtest = DATval 
    }
    # Reload best model
    load(bestmodelfile)
    #===============================================================
    # evaluate on test set
    test_factors = DATtest[,fes]
    evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
  }
}
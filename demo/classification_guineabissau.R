rm(list=ls())
graphics.off()
library(emotivepilepsy)

setwd("/home/vincent/utrecht")
shareddrive = "/media/windows-share/EEG"

funcfiles = list.files("emotivepilepsy/R",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) source(i)

trainbestmodel = TRUE #option to turn this off for Nigeria


logfile = "features_and_bestmodels/log_guinneabissau.csv" # not used when uselog = FALSE
for (aggregateperid in c(TRUE)) {
  for (proto_i in  2:1) { #"open" =1 #closed= 2
    load(file=paste0(shareddrive,"/features_and_bestmodels/features_ginneabissau_10.RData")); logdur = 10
    print("============================")
    print(paste0("protocol: ",proto_i))
    print(paste0("aggregated per id: ",aggregateperid))
    
    
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
      country = "gb" #"gb" =  guinea bisau
      bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",aggregateperid,".RData")

      save(best_model,fes,file=bestmodelfile) # Save best model
      rm(best_model,fes)
    } else {
      country = "ni" # load model from other country
      bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",logdur,"_country",country,"_perid",aggregateperid,".RData")
      LABtest = LABval # LABtest = rbind(LABval) # ignore test data, and only evaluate on training and validation data
      DATtest = DATval # DATtest = rbind(DATval)
    }
    # Reload best model
    load(bestmodelfile)
    #===============================================================
    # evaluate on test set
    test_factors = DATtest[,fes]
    print(names(DATtest[,fes])[1:10])
    
    evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
  }
}
rm(list=ls())
graphics.off()
library(emotivepilepsy)

setwd("/home/vincent/utrecht")
shareddrive = "/media/windows-share/EEG"
funcfiles = list.files("emotivepilepsy/R",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) source(i)

trainbestmodel = TRUE #option to turn this off for Nigeria
for (aggregateperid in c(TRUE)) { #FALSE
  for (proto_i in  2:1) { #"open" =1 #closed= 2
    evse = c()
    seeds2try = seq(100,1000,by=100)
    for (seed in seeds2try) { #try five seeds and select the median performing model in the test set for replication in other country
      logfile = "features_and_bestmodels/log_guinneabissau.csv" # not used when uselog = FALSE
      load(file=paste0(shareddrive,"/features_and_bestmodels/features_ginneabissau_10.RData")); logdur = 10
      print("=====")
      print(paste0("protocol: ",proto_i))
      print(paste0("seed: ",seed))
      print(paste0("aggregated per id: ",aggregateperid))
      # tidy up formatting to be suitable for classifier training
      RDL = reformat_DATLAB(DAT,LAB,aggregateperid=aggregateperid) # aggregate per unique id
      DAT =RDL$DAT
      LAB = RDL$LAB
      
      #===============================================================
      # split data in training, validation and test set
      P = split_data(LAB,DAT,logfile = logfile,proto_i=proto_i,split=c(20,20),
                     uselog = FALSE,logdur=logdur,seed=seed)
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
        bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",
                               logdur,"_country",country,"_perid",aggregateperid,"_seed",seed,".RData")
        winningmodel = trainingresults$winningmodel
        save(best_model,fes,winningmodel,file=bestmodelfile) # Save best model
      } else { # Reload best model
        country = "ni" # load model from other country
        for (k in  1: length(seeds2try)) {
          bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",
                                 logdur,"_country",country,"_perid",aggregateperid,"_seed",seeds2try[k],".RData")
          if (file.exists(bestmodelfile)) load(bestmodelfile)
        }
        LABtest = LABval # LABtest = rbind(LABval) # ignore test data, and only evaluate on training and validation data
        DATtest = DATval # DATtest = rbind(DATval)
      }
      #===============================================================
      # evaluate model on test set
      test_factors = DATtest[,fes] #DATtest[,fes]
      evaluation = evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
      evaluation$country = "gb"
      evaluation$trainingcountry = country
      evaluation$seed = seed
      evaluation$winningmodel = winningmodel
      evse = rbind(evse,evaluation)
    }
    #===============================================================
    # evaluate which seed resulted in median performance on test set and delete other models
    # we will use this model for external validation in other country
    evse = as.data.frame(evse,row.names = make.names(1:nrow(evse)))
    evse <- as.data.frame(lapply(evse, unlist))
    evse = evse[with(evse,order(test.sens,test.kappa,test.auc,test.acc)),]

    medianseed = evse$seed[round(nrow(evse)/2)]
    seeds2delete = seeds2try[which(seeds2try != medianseed)]
    for (g in 1:length(seeds2delete)) {
      file.remove(paste0(shareddrive,"/features_and_bestmodels/bestmodel_",proto_i,"_dur",
                         logdur,"_country",country,"_perid",aggregateperid,"_seed",seeds2delete[g],".RData"))
    }
    print(evse)
    # store seed evaluation
    save(evse,file=paste0(shareddrive,"/features_and_bestmodels/seedcomparison_",proto_i,"_dur",
                               logdur,"_country",country,"_perid",aggregateperid,".RData"))
  }
}

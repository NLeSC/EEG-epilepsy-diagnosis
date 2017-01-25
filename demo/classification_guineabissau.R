rm(list=ls())
graphics.off()
library(emotivepilepsy)

setwd("/home/vincent/utrecht")
shareddrive = "/media/windows-share/EEG"
funcfiles = list.files("emotivepilepsy/R",include.dirs=TRUE,full.names = TRUE)
for (i in funcfiles) source(i)

trainbestmodel = FALSE
for (epochlength in c(4)) { # in seconds
  for (aggregateperid in c(TRUE,FALSE)) { #FALSE
    for (proto_i in  2:1) { #"open" =1 #closed= 2
      evse = c()
      seeds2try = seq(100,1000,by=50)
      for (seed in seeds2try) { #try five seeds and select the median performing model in the test set for replication in other country
        load(file=paste0(shareddrive,"/features_and_bestmodels/features/features_ginneabissau_",epochlength,".RData"))
        print("=====")
        print(paste0("protocol: ",proto_i," seed: ",seed," aggregated per id: ",aggregateperid))
        # tidy up formatting to be suitable for classifier training
        RDL = reformat_DATLAB(DAT,LAB,aggregateperid=aggregateperid) # aggregate per unique id
        DAT =RDL$DAT
        LAB = RDL$LAB
        #===============================================================
        # split data in training, validation and test set
        P = split_data(LAB,DAT,proto_i=proto_i,split=c(20,20),seed=seed)
        LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
        testforoverlap = length(unique(c(LABval$id,LABtrain$id,LABtest$id))) == length(c(unique(LABval$id),unique(LABtrain$id),unique(LABtest$id)))
        if (testforoverlap == FALSE) stop("Error: Matching id numbers between subsets")
        # generate dictionary of model characteristics
        modeldict = create_modeldict(DAT)
        #===============================================================
        # train models or loaded previously trained model
        if (trainbestmodel == TRUE) {
          trainingresults = train_model(DATtrain,LABtrain,DATval,LABval,modeldict,classifier = "rf") #,aggregateperid=aggregateperid
          best_model = trainingresults$best_model
          modelcomparison = trainingresults$result
          fes = c(trainingresults$fes,which(names(DATtest) %in% c("id","diagnosis","protocol") == TRUE))
          country = "gb"
          bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodels/bestmodel_",proto_i,"_dur",
                                 epochlength,"_country",country,"_perid",aggregateperid,"_seed",seed,".RData")
          winningmodel = trainingresults$winningmodel
          save(best_model,fes,winningmodel,file=bestmodelfile) # Save best model
          #===============================================================
          # evaluate model on test set
          test_factors = DATtest[,fes]
          evaluation = evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
          evaluation$country = "gb"
          evaluation$trainingcountry = country
          evaluation$seed = seed
          evaluation$winningmodel = winningmodel
          evse = rbind(evse,evaluation)
        } else { # Reload best model
          country = "ni" # load model from other country
          for (k in  1: length(seeds2try)) {
            bestmodelfile = paste0(shareddrive,"/features_and_bestmodels/bestmodels/bestmodel_",proto_i,"_dur",
                                   epochlength,"_country",country,"_perid",aggregateperid,"_seed",seeds2try[k],".RData")
            if (file.exists(bestmodelfile)) {
              load(bestmodelfile)
              #===============================================================
              # evaluate model on test set
              test_factors = DATtest[,fes]
              evaluation = evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
              evaluation$country = "gb"
              evaluation$trainingcountry = country
              evaluation$seed = seed
              evaluation$winningmodel = winningmodel
              evse = rbind(evse,evaluation)
            }
          }
        }
        # #===============================================================
        # # evaluate model on test set
        # test_factors = DATtest[,fes]
        # evaluation = evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
        # evaluation$country = "gb"
        # evaluation$trainingcountry = country
        # evaluation$seed = seed
        # evaluation$winningmodel = winningmodel
        # evse = rbind(evse,evaluation)
      }
      #===============================================================
      # evaluate which seed resulted in median performance on test set and delete other models
      # we will use this model for external validation in other country
      evse = as.data.frame(evse,row.names = make.names(1:nrow(evse)))
      evse <- as.data.frame(lapply(evse, unlist))
      evse = evse[with(evse,order(test.acc,test.kappa,test.auc,test.sens)),] #val.acc,val.kappa,val.auc,val.sens
      print(evse)
      # store seed evaluation
      save(evse,file=paste0(shareddrive,"/features_and_bestmodels/evaluation/seedcomparison_",proto_i,"_dur",
                            epochlength,"_country",country,"_perid",aggregateperid,"_evalcountry",evaluation$country,".RData"))
    }
  }
}
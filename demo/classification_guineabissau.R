rm(list=ls())
#==============================================
# Update the following lines:
setwd("/home/vincent/utrecht")
outputdir = "/media/windows-share/EEG" # this is where folders will be created to store the output
namecountry = "gb"
namecountry2 = "ni" #the country for cross validation

#==============================================
funcfiles = list.files("EEG-epilepsy-diagnosis/R",include.dirs=TRUE,full.names = TRUE) # this line only needed when developing
for (i in funcfiles) source(i) # this line only needed when developing

outputdir_features = paste0(outputdir,"/EEGs_features") # directiory should have been generated in the pre-processing
outputdir_bestmodels = paste0(outputdir,"/EEGs_bestmodels")
if (!file.exists(outputdir_bestmodels)) dir.create(outputdir_bestmodels)
outputdir_evaluation =  paste0(outputdir,"/EEGs_evaluation")
if (!file.exists(outputdir_evaluation)) dir.create(outputdir_evaluation)

trainbestmodel = FALSE
limit2sdfeatutes = TRUE
for (epochlength in c(4)) { # in seconds
  for (aggregateperid in c(FALSE,TRUE)) { #FALSE
    for (proto_i in  2) { #"open" =1 #closed= 2
      evse = c()
      seeds2try = seq(100,600,by=50)
      for (seed in seeds2try) { #try five seeds and select the median performing model in the test set for replication in other country
        load(file=paste0(outputdir_features,"/features_",namecountry,"_epoch",epochlength,".RData"))
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
        featuredict = create_featuredict(DAT)
        
        if (limit2sdfeatutes == TRUE) {
          # ok, let's see how well the model performans if we only keep the features related to sd which have been showing up to dominate model performance
          featurenamestokeep = rownames(featuredict[which(featuredict$feature == "sd"),])
          col2keep = which(colnames(DAT) %in% c(featurenamestokeep,"id","diagnosis","protocol") == TRUE)
          DATtrain = DATtrain[,col2keep]
          DATval = DATval[,col2keep]
          DATtest = DATtest[,col2keep]
          featuredict = featuredict[which(featuredict$feature == "sd"),]
        }
        kkk
        #===============================================================
        # train models or loaded previously trained model
        if (trainbestmodel == TRUE) {
          trainingresults = train_model(DATtrain,LABtrain,DATval,LABval,featuredict,classifier = "rf") #,aggregateperid=aggregateperid
          best_model = trainingresults$best_model
          modelcomparison = trainingresults$result
          fes = c(trainingresults$fes,which(names(DATtest) %in% c("id","diagnosis","protocol") == TRUE))
          bestmodelfile = paste0(outputdir_bestmodels,"/bestmodel_",proto_i,"_dur",
                                 epochlength,"_country",namecountry,
                                 "_perid",aggregateperid,"_SDonly",limit2sdfeatutes ,
                                 "_seed",seed,".RData")
          winningmodel = trainingresults$winningmodel
          save(best_model,fes,winningmodel,file=bestmodelfile) # Save best model
          #===============================================================
          # evaluate model on test set
          test_factors = DATtest[,fes]
          evaluation = evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
          evaluation$country = namecountry
          evaluation$trainingcountry = namecountry
          evaluation$seed = seed
          evaluation$winningmodel = winningmodel
          evse = rbind(evse,evaluation)
        } else { # Reload best model from the other country
          for (k in  1: length(seeds2try)) {
            bestmodelfile = paste0(outputdir_bestmodels,"/bestmodel_",proto_i,"_dur",
                                   epochlength,"_country",namecountry2,"_perid",aggregateperid,
                                   "_SDonly",limit2sdfeatutes ,"_seed",seeds2try[k],".RData")
            if (file.exists(bestmodelfile)) {
              load(bestmodelfile)
              #===============================================================
              # evaluate model on test set
              test_factors = DATtest[,fes]
              evaluation = evaluatemodel(model=best_model,x=test_factors,labels=LABtest,proto_i=proto_i,aggregateperid=aggregateperid)
              evaluation$country = namecountry
              evaluation$trainingcountry = namecountry2
              evaluation$seed = seed
              evaluation$winningmodel = winningmodel
              evse = rbind(evse,evaluation)
            }
          }
        }
      }
      #===============================================================
      # Order and store performance of 11 models across 11 seeds (if done in other country then this results in 11 x 11 evaluations)
      evse = as.data.frame(evse,row.names = make.names(1:nrow(evse)))
      evse <- as.data.frame(lapply(evse, unlist))
      evse = evse[with(evse,order(test.acc,test.kappa,test.auc,test.sens)),] #val.acc,val.kappa,val.auc,val.sens
      save(evse,file=paste0(outputdir_evaluation,"/seedcomparison_",proto_i,"_dur",
                            epochlength,"_country",evaluation$trainingcountry,"_perid",aggregateperid,
                            "_evalcountry",evaluation$country,"_SDonly",limit2sdfeatutes,".RData"))
    }
  }
}
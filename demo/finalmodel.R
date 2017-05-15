rm(list=ls())
#============================================
# Script to train final model on all data from GB and evaluate it in all data from NI
#==============================================
# Update the following lines:
setwd("/home/vincent/utrecht")
outputdir = "/media/windows-share/EEG" # this is where folders will be created to store the output
namecountry = "gb"
namecountry2 = "ni" #the country for cross validation


# What was most popular wavelet?
load(paste0(outputdir,"/EEGs_evaluation/seedcomparison_2_dur4_countrygb_peridTRUE_evalcountrygb_SDonlyTRUE.RData"))
print("Winnning models across seeds")
winners = table(evse$winningmodel)
print(winners)
bestwavelet = names(winners[which.max(winners)[1]])

#==============================================
# get functions and folders:
funcfiles = list.files("EEG-epilepsy-diagnosis/R",include.dirs=TRUE,full.names = TRUE) # this line only needed when developing
for (i in funcfiles) source(i) # this line only needed when developing
outputdir_features = paste0(outputdir,"/EEGs_features") # directiory should have been generated in the pre-processing
outputdir_finalmodel = paste0(outputdir,"/EEGs_finalmodel")
if (!file.exists(outputdir_finalmodel)) dir.create(outputdir_finalmodel)
outputdir_finalmodelevaluation =  paste0(outputdir,"/EEGs_finalmodelevaluation")
if (!file.exists(outputdir_finalmodelevaluation)) dir.create(outputdir_finalmodelevaluation)
#==============================================
trainbestmodel = TRUE
limit2sdfeatutes = TRUE
proto_i = 2
aggregateperid = TRUE
epochlength = 4
performancemetric = "Accuracy"
TL = 3
# Load GB data
load(file=paste0(outputdir_features,"/features_",namecountry,"_epoch",epochlength,".RData"))
# tidy up formatting to be suitable for classifier training
RDL = reformat_DATLAB(DAT,LAB,aggregateperid=aggregateperid) # aggregate per unique id
DAT =RDL$DAT
LAB = RDL$LAB
featuredict = create_featuredict(DAT)
if (limit2sdfeatutes == TRUE) {
  # ok, let's see how well the model performans if we only keep the features related to sd which have been showing up to dominate model performance
  featurenamestokeep = rownames(featuredict[which(featuredict$feature == "sd"),])
  col2keep = which(colnames(DAT) %in% c(featurenamestokeep,"id","diagnosis","protocol") == TRUE)
  DAT = DAT[,col2keep]
  featuredict = featuredict[which(featuredict$feature == "sd"),]
}
# Train models
print("Train model")
ctrl = caret::trainControl(method = "repeatedcv",number=5,repeats=1,search="random",
                           classProbs = TRUE) #, summaryFunction = caret::prSummary) #twoClassSummary)  # twoClassSummary is needed for Sensitivity optimization
allvalues = featuredict$wvtype
valueseval = unique(allvalues)
fes = which(allvalues==valueseval[which(valueseval == bestwavelet)] | allvalues == "raw") # always add features derived from raw data
train_factors = DAT[,fes]
m_rf = caret::train(y=as.factor(make.names(DAT$diagnosis)),x=train_factors,seeds=100,
                    method="rf",importance=TRUE,trControl=ctrl,tuneLength=TL,metric=performancemetric) # # train 5 different mtry values using random search

rm(LAB, DAT)
#==============================================
# Get Nigeria data:
#==============================================
print("Evaluate model in Nigerian data")
load(file=paste0(outputdir_features,"/features_",namecountry2,"_epoch",epochlength,".RData"))
LAB$diagnosis = as.character(LAB$diagnosis)
LAB$diagnosis[which(LAB$diagnosis == "control")] = "Control"
LAB$diagnosis[which(LAB$diagnosis == "epilepsy")] = "Epilepsy"
LAB$diagnosis = as.factor(LAB$diagnosis)
# tidy up formatting to be suitable for classifier training
RDL = reformat_DATLAB(DAT,LAB,aggregateperid=aggregateperid) # aggregate per unique id
DAT =RDL$DAT
LAB = RDL$LAB
featuredict = create_featuredict(DAT)
if (limit2sdfeatutes == TRUE) {
  # ok, let's see how well the model performans if we only keep the features related to sd which have been showing up to dominate model performance
  featurenamestokeep = rownames(featuredict[which(featuredict$feature == "sd"),])
  col2keep = which(colnames(DAT) %in% c(featurenamestokeep,"id","diagnosis","protocol") == TRUE)
  DAT = DAT[,col2keep]
  featuredict = featuredict[which(featuredict$feature == "sd"),]
}
fes = c(fes,which(names(DAT) %in% c("id","diagnosis","protocol") == TRUE))

# apply model:
pred_test = predict(m_rf,DAT[,fes],type="prob")
evaluation = evaluatemodel(model=m_rf,x=DAT[,fes],labels=LAB,proto_i=proto_i,aggregateperid=aggregateperid)



# Create roc curve
labels_agg = LAB
x_agg = DAT[,fes]
# result.roc <- pROC::roc(x_agg$diagnosis, pred_test$X1)
# auctest = result.roc$auc
# result.coords <- pROC::coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
# pred_test_cat = rep("X1",nrow(pred_test))
# pred_test_cat[pred_test$X2 > 0.500] = "X2"
# refe = make.names(labels_agg$diagnosis)

library(pROC)
rocobj = pROC::roc(x_agg$diagnosis-1, pred_test$X2, direction="<")
ci.sp.obj <- ci.sp(rocobj, sensitivities=seq(0, 1, .01), boot.n=200)

outfile = paste0(outputdir,paste0("/Figure2_ROC.jpeg"))
jpeg(filename=outfile, units="in",width = 3.5,height= 3.5,res=600,pointsize = 12)
par(bty="l",cex.axis=0.9,cex=0.9)

plot(rocobj,print.auc=TRUE,lwd=3, bty="l") # restart a new plot
plot(ci.sp.obj, type="shape", col="gray")

dev.off()


# Look for associations with descriptives:
NIMETA = read.csv("/media/windows-share/EEG/input/merged-meta-data_nigeria.csv")
N = data.frame(Pred=pred_test$X2,id=DAT$id,Truth=x_agg$diagnosis-1) # 2 is Epilepsy
N2 = merge(N,NIMETA,by.x="id",by.y="subject.id")

NDQ = read.csv("/media/windows-share/EEG/quality_overview_EEGs_ni_cleaned.csv")
cut = which(is.na(NDQ$id) == TRUE)
if (length(cut) > 0) NDQ = NDQ[-cut,]
# normalize scores
NDQ$totalqc = NDQ$closed_qc1 + NDQ$closed_qc2 + NDQ$closed_qc3 + NDQ$closed_qc4
NDQ$closed_qc1 = NDQ$closed_qc1 / NDQ$totalqc
NDQ$closed_qc2 = NDQ$closed_qc2 / NDQ$totalqc
NDQ$closed_qc3 = NDQ$closed_qc3 / NDQ$totalqc
NDQ$closed_qc4 = NDQ$closed_qc4 / NDQ$totalqc

NDQ$totalqu = NDQ$closed_quality0 + NDQ$closed_quality1 + NDQ$closed_quality2
NDQ$closed_quality0 = NDQ$closed_quality0 / NDQ$totalqu
NDQ$closed_quality1 = NDQ$closed_quality1 / NDQ$totalqu
NDQ$closed_quality2 = NDQ$closed_quality2 / NDQ$totalqu

N3 = merge(N2,NDQ,by="id")


N3$succes = abs(round(N3$Pred)-N3$Truth)
fit = glm(succes ~ age+ sex+closed_mean_RAWCQ+closed_quality2+closed_quality1,family=binomial(link='logit'),data=N3)
print(summary(fit))





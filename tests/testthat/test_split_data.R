library(emotivepilepsy)
test_that("split_data splits dataset in training validation and test set", {
  data(data.eeg)
  data(data.labels)
  RDL = reformat_DATLAB(data.eeg,data.labels,aggregateperid=TRUE)
  DAT =RDL$DAT
  LAB = RDL$LAB
  P = split_data(LAB,DAT,logfile = logfile,proto_i=1,split=c(2,2),
                 uselog = FALSE,logdur=10)
  LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  expect_that(LABval,is_a("data.frame"))
})
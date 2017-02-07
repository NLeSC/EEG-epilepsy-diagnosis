library(emotivepilepsy)
test_that("create_featuredict creates dictionary of possible models to consider", {
  data(data.eeg)
  data(data.labels)
  aggregateperid = TRUE
  RDL = reformat_DATLAB(data.eeg,data.labels,aggregateperid=aggregateperid)
  DAT =RDL$DAT
  LAB = RDL$LAB
  P = split_data(LAB,DAT,proto_i=1,split=c(1,1),seed=300)
  LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  featuredict = create_featuredict(DAT)
  expect_that(featuredict,is_a("data.frame"))
})
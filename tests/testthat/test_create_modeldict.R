library(emotivepilepsy)
test_that("create_modeldict creates dictionary of possible models to consider", {
  data(data.eeg)
  data(data.labels)
  aggregateperid = TRUE
  RDL = reformat_DATLAB(data.eeg,data.labels,aggregateperid=aggregateperid)
  DAT =RDL$DAT
  LAB = RDL$LAB
  P = split_data(LAB,DAT,proto_i=1,split=c(1,1),seed=300)
  LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  modeldict = create_modeldict(DAT)
  expect_that(modeldict,is_a("data.frame"))
})
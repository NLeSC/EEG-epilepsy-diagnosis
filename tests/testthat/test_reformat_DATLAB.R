library(emotivepilepsy)
test_that("reformat_DATLAB reformats the dataframes eegdata and eegdatalabels", {
  data(data.eeg)
  data(data.labels)
  RDL = reformat_DATLAB(data.eeg,data.labels,aggregateperid=TRUE)
  expect_that(RDL$DAT,is_a("data.frame"))
})
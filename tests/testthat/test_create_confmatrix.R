library(emotivepilepsy)
test_that("create_confmatrix produces a table", {
  set.seed(21)
  reference <- make.names(round(runif(n=10,min=0,max=0.8),digits=0))
  set.seed(23)
  predicted <- make.names(round(runif(n=10,min=0,max=0.8),digits=0))
  confusionmatrix <- create_confmatrix(reference, predicted)
  expect_that(confusionmatrix,is_a("table"))
})
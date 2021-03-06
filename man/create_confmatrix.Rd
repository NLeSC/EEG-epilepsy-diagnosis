% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_confmatrix.R
\name{create_confmatrix}
\alias{create_confmatrix}
\title{Create confusion matrix}
\usage{
create_confmatrix(predicted, reference)
}
\arguments{
\item{predicted}{predicted class}

\item{reference}{reference class}
}
\value{
confmat confusion matrix
}
\description{
Creates an NxN confusion matrix for a binary classification and its reference values
}
\examples{
set.seed(21)
reference <- make.names(round(runif(n=10,min=0,max=0.8),digits=0))
set.seed(23)
predicted <- make.names(round(runif(n=10,min=0,max=0.8),digits=0))
confusionmatrix <- create_confmatrix(reference, predicted)
}
\author{
Vincent T van Hees
}

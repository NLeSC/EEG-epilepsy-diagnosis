\name{evaluatemodel}
\alias{evaluatemodel}
\title{Evaluate model on data}
\usage{
evaluatemodel(model,x,labels,proto_i,aggregateperid)
}
\arguments{
\item{model}{Model produced by caret function train}

\item{x}{Test data on which the model will be applied}

\item{labels}{Data labels}

\item{proto_i}{Protocol, 1 = eyes open, 2 = eyes closed, this is only used
for logging and not for the evaluation}

\item{aggregateperid}{Boolean to indicate whether features were aggregated
per participant (TRUE) or not (FALSE)}
}

\value{
\item{\code{proto_i}}{same as input}
\item{\code{test.acc}}{Test set accuracy}
\item{\code{test.kappa}}{Test set Cohen's kappa coeficient}
\item{\code{test.auc}}{Test set Area under the curve}
\item{\code{test.confmatrix}}{Test set confusion matrix}
\item{\code{test.sens}}{Test set sensitivity}
\item{\code{aggregateperid}}{same as input}
}
\description{
Evaluates a model on test dataset
}
\author{
Vincent T van Hees
}

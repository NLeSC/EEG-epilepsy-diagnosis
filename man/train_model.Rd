\name{train_model}
\alias{train_model}
\title{Train binary classifier}
\usage{
train_model(DATtrain,LABtrain,DATval,LABval,modeldict,classifier="rf")
}
\arguments{
\item{DATtrain}{Training feature values as produced by extract_features and split 
up by split_data}

\item{LABtrain}{Training labels as derived from split_data}

\item{DATval}{Validation feature values as produced by extract_features and split
up by split_data}

\item{LABval}{Validation labels as derived from split_data}

\item{modeldict}{Dataframe with all overview of all available combinations of wavelet
types, wavelet levels, extracted signal features, and aggregation type.}

\item{classifier}{type of classification model, defaul is "rf" for random forrest.
Alternative value could be "lg" for logistic regression}


\value{
\item{\code{result}}{classificiation performance in the validation set for each value
of the hyperparameter}
\item{\code{best_model}}{classification model as produced by caret function train}
\item{\code{fes}}{column indices of the DAT corresponding to the features used for the best model}
}
\description{
Loads EEG data and writes cleaned data to an output directory
}
\examples{
\dontrun{
  data(data.eeg)
  data(data.labels)
  RDL = reformat_DATLAB(data.eeg,data.labels,aggregateperid=TRUE)
  DAT =RDL$DAT
  LAB = RDL$LAB
  P = split_data(LAB,DAT,logfile = logfile,proto_i=1,split=c(2,2),
                 uselog = FALSE,logdur=10)
  LABval=P$LABval;LABtest=P$LABtest;LABtrain=P$LABtrain;DATval=P$DATval;DATtest=P$DATtest;DATtrain=P$DATtrain
  modeldict = create_modeldict(DAT)
  mymodel = train_model(DATtrain,LABtrain,DATval,LABval,modeldict,classifier="rf")
}
}
\author{
Vincent T van Hees
}

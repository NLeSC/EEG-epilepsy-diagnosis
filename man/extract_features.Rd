\name{extract_features}
\alias{extract_features}
\title{Creates a dictionary (dataframe) of all features in the data}
\usage{
extract_features(datadir,sf,n.levels,filtertypes,epochlength,fn)
}
\arguments{
\item{datadir}{location of the csv-files with cleaned data per person}
\item{sf}{sample frequency in Hertz}
\item{n.levels}{number of wavelet levels to extract}
\item{filtertypes}{wavelet filter types as expected by the mra function from the wavelets package, e.g. c("d2","d4")}
\item{epochlength}{duration of epoch (time window) in seconds}
\item{fn}{names of features to extract from wavelets, e.g. "sd" to extract the standardard. Underlying code
includes explorative work on other features}
}

\value{
\item{DAT}{dataframe with features}
\item{LAB}{dataframe with labels}
}
\description{
Extracts features from raw EEG data using wavelets
}
\author{
Vincent T van Hees
}

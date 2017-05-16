\name{clean_emotiv}
\alias{clean_emotiv}
\title{Clean emotive EEG data}
\usage{
clean_emotiv(datadir, metadatafile,  outputdir, knownerrors, protocoltimes,
referencegroup, condition_start_closed, protocolvariable, outputdir_clean,
sf=128, gyrothreshold=30, mindur=4)
}
\arguments{
\item{datadir}{data directory}

\item{metadatafile}{file with meta}

\item{outputdir}{directory where output will be stored}

\item{knownerrors}{data frame with 3 columns for the person identifier,
start time and end time corresponding to the data that needs to be ignored
Time format is minute:seconds}

\item{protocoltimes}{vector of length 3 with start times of first and 
  second experiment and end time of last experiment}

\item{condition_start_closed}{protocol label corresponding to the protocol
  where the start was with eyes closed, e.g. "closed-3min-then-open-2min"}

\item{referencegroup}{diagnosis label corresponding to the reference group, e.g. "Control"}

\item{protocolvariable}{ Name of column in metadata that indicates the 
experimental protocol, e.g. "Eyes.condition"}

\item{outputdir_clean}{ Directory for storage of cleaned data}

\item{sf}{sample frequency}

\item{gyrothreshold}{threshold for gyroscope signal above which data is considered not clean}

\item{mindur}{required minimum duration (minutes) of the EEG time series}


}
\value{
matrix with overview of how much data was available per participant id
}
\description{
Loads EEG data and writes cleaned data to an output directory
}
\examples{
\dontrun{
datafile = "C:/myfolder/example.csv"
clean_emotiv(datadir = "C:/mydatafolder", metadatafile="C:/metadata.csv",
  outputdir="C:/output", sf=128, gyrothreshold=30, mindur=4,
  knownerrors=list(1,"1:10","1:20"), protocoltimes=c(60,180,300),
  referencegroup="Control", condition_start_closed="closed-3min-then-open-2min")
}
}
\author{
Vincent T van Hees
}


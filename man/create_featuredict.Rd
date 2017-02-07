\name{create_featuredict}
\alias{create_featuredict}
\title{Creates a dictionary (dataframe) of all features in the data}
\usage{
create_featuredict(DAT)
}
\arguments{
\item{DAT}{Dataframe DAT produced by extract_features}
}

\value{
\item{\code{featuredict}}{dataframe with columns for waveletype, feature, waveletlevel and aggregation per id indicator}
}
\description{
Creates a dataframe of all features in the data with one characteristic per column
}
\author{
Vincent T van Hees
}

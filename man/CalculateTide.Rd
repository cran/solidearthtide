\name{CalculateTide}
\alias{CalculateTide}
\title{
Predicted solid earth tide displacement
}
\description{
\code{CalculateTide} Predicteds the solid earth tide displacement.
}
\usage{
CalculateTide(fname, yearIn, monthIn, dayIn, latIn, lonIn, boolRet=FALSE)
}
\arguments{
  \item{fname}{ optional filename for output}
  \item{yearIn}{year}
  \item{monthIn}{month}
  \item{dayIn}{day}
  \item{latIn}{latitude}
  \item{lonIn}{longitude}
  \item{boolRet}{ optional boolean, TRUE=return data as a matrix}
}
\examples{
latIn <- 30.25
lonIn <- -97.75
yearIn <- 2015
monthIn <- 5
dayIn <- 11
CalculateTide('',yearIn, monthIn, dayIn, latIn, lonIn, TRUE)
}
\author{
Jose Gama
}
\keyword{ programming }

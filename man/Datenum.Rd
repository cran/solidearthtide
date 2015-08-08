\name{Datenum}
\alias{Datenum}
\title{
Convert date to serial date number as Matlab
}
\description{
\code{Datenum} converts date to serial date number with 
the same values as Matlab's datenum
}
\usage{
Datenum(d)
}
\arguments{
  \item{d}{ date object}
}
\value{serial date number
}
\examples{
Datenum(as.Date("may 27, 1984", "\%b \%d, \%Y"))
}
\author{
Jose Gama
}
\keyword{ programming }

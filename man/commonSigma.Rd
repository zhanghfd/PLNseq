\name{commonSigma}
\alias{commonSigma}

\title{
Common sigma
}
\description{
Estimate `mu' (mean parameter lognormal distribution) for each gene and condition and a common `sigma' (standard deviation parameter of lognormal distribution).
}
\usage{
commonSigma(d)
}

\arguments{
  \item{d}{
  This is a PLNseq object.
  }
}

\value{
    \item{d$commonSigma}{A common `sigma'}
}

\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,2);
    d = sizeFactor(d,maxCount=2e3);
    d = commonSigma(d);

}


\keyword{mean}
\keyword{standard deviation}
\keyword{lognormal distribution}

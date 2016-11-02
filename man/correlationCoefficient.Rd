\name{correlationCoefficient}
\alias{correlationCoefficient}

\title{
Correlation coefficient
}
\description{
Estimate correlation coefficient parameter and its standard error in the multivariate lognormal distribution.
}

\usage{
correlationCoefficient(d)
}

\arguments{
  \item{d}{
  This is a PLNseq object.
  }
}

\value{
    \item{d$rho}{Correlation coefficient `rho' in the multivariate lognormal distribution}
    \item{d$rho.se}{Standard error of estimated `rho'}
}

\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,2);
    d = sizeFactor(d,maxCount=2e3);
    d = commonSigma(d);
    d = correlationCoefficient(d);
}

\keyword{correlation coefficient} 
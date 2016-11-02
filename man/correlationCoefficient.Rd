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
    
    ## common correlation
	## commonCorrelation = TRUE;
    ## d = correlationCoefficient(d);
    
    ## clustered correlation
	## commonCorrelation = FALSE;
    ## J = nrow(count);
    ## J1 = round(J/2);
    ## d$cluster = c(rep(1,J1),rep(2,J-J1));
    ## d = correlationCoefficient(d);
    
}

\keyword{correlation coefficient} 
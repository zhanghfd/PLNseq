\name{LRtest2}
\alias{LRtest2}

\title{
Likelihood ratio test for differential expression analysis with cluster-specific correlations.
}
\description{
This function calculates log-fold changes, likelihood ratio test statistics, and p-values for a list of genes. This function should be called after cluster-specific correlations are returned by `correlationCoefficient'.
}
\usage{
    LRtest2(d,M,use.commonSigma,id)
}

\arguments{
    \item{d}{
    This is a PLNseq object.
    }
    \item{M}{
	The number of simulations used in Monte-Carlo method for calculating likelihood ratio test statistics.
    }
    \item{use.commonSigma}{
    Use common `sigma' (TRUE) or genewise `sigma' (FALSE), with default value `FALSE'.
    }
    \item{id}{
    A vector consisting of a subset of {1,...,J}, with default value 1:J.
    }
}

\value{
    \item{LR}{Estimation and test results: `log-FC', `LR statistic', `p value'.}
}



\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,conditionNumber=2);
    d = sizeFactor(d,maxCount=2e3);
    
    ##  Not run:
    ##  d = commonSigma(d);
    ##  J = nrow(count);
    ##  J1 = round(J/2);
    ##  d$commonCorrelation = FALSE;
    ##  d$cluster = c(rep(1,J1),rep(2,J-J1));
    ##  d = correlationCoefficient(d);
    ##  d = genewiseSigma(d);
    ##  d = LRtest2(d,M=3e4,use.commonSigma=FALSE,id=1:100);
}

\keyword{differential expression}
\keyword{likelihood ratio test}
\keyword{correlated read count}
\keyword{cluster-specific correlation}

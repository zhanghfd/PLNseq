\name{LRtest1}
\alias{LRtest1}

\title{
Likelihood ratio test for differential expression analysis with common correlation.
}
\description{
This function calculates log-fold changes, likelihood ratio test statistics, and p-values for a list of genes. This function should be called after a commom correlation matrix is returned by `correlationCoefficient'.
}
\usage{
    LRtest1(d,z,use.commonSigma,id)
}

\arguments{
    \item{d}{
    This is a PLNseq object.
    }
    \item{z}{
    J independent samples (a matrix of dimension J by R) drown from multivariate normal distribution with expectations 0, variances 1, and a common correlation coefficient matrix estimated by `correlationCoefficient'.
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
    ##  d$commonCorrelation = TRUE;
    ##  d = correlationCoefficient(d);
    ##  d = genewiseSigma(d);
    ##  library('MASS');
    ##  z = mvrnorm(n=1e5,mu=rep(0,2),Sigma=d$rho);
    ##  d = LRtest1(d,z,use.commonSigma=FALSE,id=1:100);
}

\keyword{differential expression}
\keyword{likelihood ratio test}
\keyword{correlated read count}
\keyword{common correlation}

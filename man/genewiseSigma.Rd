\name{genewiseSigma}
\alias{genewiseSigma}

\title{
Genewise sigma
}
\description{
Estimate genewise `sigma' (standard deviation parameter of lognormal distribution).
}
\usage{
genewiseSigma(d,w=25)
}

\arguments{
  \item{d}{
  This is a PLNseq object.
  }
  \item{w}{
  Shrinkage parameter.
  }
}

\value{
    \item{d$genewiseSigma}{Genewise `sigma'}
}

\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,2);
    d = sizeFactor(d,maxCount=2e3);

    d = genewiseSigma(d);

}


\keyword{mean}
\keyword{standard deviation}
\keyword{lognormal distribution}

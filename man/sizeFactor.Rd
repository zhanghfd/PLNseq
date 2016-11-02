\name{sizeFactor}
\alias{sizeFactor}

\title{
Estimate size factor for each sample.
}
\description{
Estimate size factor for each sample using median normalization method.
}

\usage{
sizeFactor(d,maxCount)
}

\arguments{
  \item{d}{
  This is a PLNseq object.
  }
  \item{maxCount}{
  The maximal count after shrinkage, with a default value NA (no shrinkage).
  }
}

\value{
    \item{d$sample$sizeFactor}{Estimated size factors}
}

\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,conditionNumber=2);
    d = sizeFactor(d,maxCount=2e3);

}

\keyword{median normalization}
\keyword{size factor}

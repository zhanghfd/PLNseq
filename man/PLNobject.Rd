\name{PLNobject}
\alias{PLNobject}

\title{
PLN object
}
\description{
Create a PLN object, a list containing a read count matrix `count' and sample description matrix `sample'.
}
\usage{
PLNobject(count,conditionNumber)
}

\arguments{
  \item{count}{
  This is a matrix containing the read counts of R*I samples at J genes (R is the number of conditions in each matched sample set and I is the number of sample sets). Here columns 1 through I are for I     
  independent samples from condition 1, columns I+1 through 2I are for I samples from condition 2 matched 
  by samples 1 through I, ... , columns (R-1)*I+1 through R*I are for I samples from condition R matched 
  by samples 1 through I.
  }
  \item{conditionNumber}{
  Number of conditions.
  }
}

\value{
    \item{d$count}{Origninal read count matrix}
    \item{d$conditionNumber}{The number of conditions}
    \item{d$sample}{A matrix of sample information: `SampleName', `TotalCount', `MedianCount'}
}

\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,2);

}


\keyword{RNA-seq read count}
\keyword{matched sample}
